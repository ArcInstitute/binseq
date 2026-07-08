//! `zstd-rs` backend: a thin safe wrapper around `libzstd-rs-sys`.
//!
//! `libzstd-rs-sys` is a c2rust transpilation of the zstd C library. It exposes the
//! same flat, C-ABI-shaped `ZSTD_*` functions as `zstd-sys` (same names, same `size_t`
//! error-code convention), but implemented in pure Rust with no C compiler
//! dependency. There is no safe Rust wrapper published yet (its README describes the
//! primary use today as "a static library linked into C projects"), so this module is
//! one: it mirrors the exact subset of `zstd`/`zstd_safe`'s API that the rest of this
//! crate relies on, so call sites are identical across both backends.
//!
//! All `unsafe` FFI calls are confined to this module.

use std::ffi::CStr;
use std::io::{self, Read, Write};
use std::marker::PhantomData;
use std::os::raw::c_void;

use libzstd_rs_sys as ffi;

pub(crate) mod zstd_safe {
    use super::{CStr, PhantomData, c_void, ffi};

    /// Destination buffer usable by [`CCtx::compress2`] / [`DCtx::decompress`].
    ///
    /// Mirrors `zstd_safe::WriteBuf`: implemented for `[u8]` (fixed capacity) and
    /// `Vec<u8>` (capacity-bounded write, `set_len` to the actual size on success) so
    /// call sites written against either backend behave identically.
    ///
    /// # Safety
    /// Implementors must ensure `as_mut_ptr` returns a pointer valid for
    /// `capacity()` writable bytes, and that `set_len(n)` is only a no-op-safe
    /// bookkeeping update when the first `n` bytes have actually been initialized.
    pub(crate) unsafe trait WriteBuf {
        fn capacity(&self) -> usize;
        fn as_mut_ptr(&mut self) -> *mut u8;
        /// # Safety
        /// The first `len` bytes of the buffer must have been initialized.
        unsafe fn set_len(&mut self, len: usize);
    }

    unsafe impl WriteBuf for [u8] {
        fn capacity(&self) -> usize {
            self.len()
        }
        fn as_mut_ptr(&mut self) -> *mut u8 {
            <[u8]>::as_mut_ptr(self)
        }
        unsafe fn set_len(&mut self, _len: usize) {}
    }

    unsafe impl WriteBuf for Vec<u8> {
        fn capacity(&self) -> usize {
            Vec::capacity(self)
        }
        fn as_mut_ptr(&mut self) -> *mut u8 {
            Vec::as_mut_ptr(self)
        }
        unsafe fn set_len(&mut self, len: usize) {
            unsafe { Vec::set_len(self, len) };
        }
    }

    /// Checks a raw zstd return code, converting it to `Result` the same way
    /// `zstd_safe` does (the error variant is the raw code, decoded on demand via
    /// [`get_error_name`]).
    fn check(code: usize) -> Result<usize, usize> {
        if ffi::ZSTD_isError(code) != 0 {
            Err(code)
        } else {
            Ok(code)
        }
    }

    /// Compression parameters - only the two variants this crate actually sets.
    #[derive(Clone, Copy)]
    pub(crate) enum CParameter {
        CompressionLevel(i32),
        EnableLongDistanceMatching(bool),
    }

    pub struct CCtx<'a> {
        raw: *mut ffi::ZSTD_CCtx,
        _marker: PhantomData<&'a ()>,
    }

    // SAFETY: `raw` is an owned, exclusively-held heap allocation created by
    // `ZSTD_createCCtx`. It is never shared across threads without synchronization
    // in this codebase (each `CCtx` is only ever accessed through `&mut self`),
    // matching the same unsafe impls `zstd_safe::CCtx` provides.
    unsafe impl Send for CCtx<'_> {}
    unsafe impl Sync for CCtx<'_> {}

    impl CCtx<'_> {
        #[must_use]
        pub(crate) fn create() -> Self {
            let raw = unsafe { ffi::ZSTD_createCCtx() };
            assert!(!raw.is_null(), "ZSTD_createCCtx returned null");
            Self {
                raw,
                _marker: PhantomData,
            }
        }

        pub(crate) fn set_parameter(&mut self, param: CParameter) -> Result<usize, usize> {
            let (param, value) = match param {
                CParameter::CompressionLevel(level) => {
                    (ffi::ZSTD_cParameter::ZSTD_c_compressionLevel, level)
                }
                CParameter::EnableLongDistanceMatching(enable) => (
                    ffi::ZSTD_cParameter::ZSTD_c_enableLongDistanceMatching,
                    i32::from(enable),
                ),
            };
            check(unsafe { ffi::ZSTD_CCtx_setParameter(self.raw, param, value) })
        }

        pub(crate) fn compress2<C: WriteBuf + ?Sized>(
            &mut self,
            dst: &mut C,
            src: &[u8],
        ) -> Result<usize, usize> {
            let written = check(unsafe {
                ffi::ZSTD_compress2(
                    self.raw,
                    dst.as_mut_ptr().cast::<c_void>(),
                    dst.capacity(),
                    src.as_ptr().cast::<c_void>(),
                    src.len(),
                )
            })?;
            unsafe { dst.set_len(written) };
            Ok(written)
        }
    }

    impl Drop for CCtx<'_> {
        fn drop(&mut self) {
            unsafe { ffi::ZSTD_freeCCtx(self.raw) };
        }
    }

    pub struct DCtx<'a> {
        raw: *mut ffi::ZSTD_DCtx,
        _marker: PhantomData<&'a ()>,
    }

    // SAFETY: see the matching comment on `CCtx`'s `Send`/`Sync` impls above.
    unsafe impl Send for DCtx<'_> {}
    unsafe impl Sync for DCtx<'_> {}

    impl DCtx<'_> {
        #[must_use]
        pub(crate) fn create() -> Self {
            let raw = unsafe { ffi::ZSTD_createDCtx() };
            assert!(!raw.is_null(), "ZSTD_createDCtx returned null");
            Self {
                raw,
                _marker: PhantomData,
            }
        }

        pub(crate) fn decompress<C: WriteBuf + ?Sized>(
            &mut self,
            dst: &mut C,
            src: &[u8],
        ) -> Result<usize, usize> {
            let written = check(unsafe {
                ffi::ZSTD_decompressDCtx(
                    self.raw,
                    dst.as_mut_ptr().cast::<c_void>(),
                    dst.capacity(),
                    src.as_ptr().cast::<c_void>(),
                    src.len(),
                )
            })?;
            unsafe { dst.set_len(written) };
            Ok(written)
        }

        /// Wraps `ZSTD_decompressStream`. Returns `(bytes written to dst, bytes
        /// consumed from src, hint)`, where `hint == 0` means the current frame is
        /// complete (matching the C API's contract). Used by [`super::stream::copy_decode`]
        /// so decoding doesn't depend on the frame carrying an embedded content
        /// size (streaming-produced frames - including ones from older files on
        /// disk - often don't).
        pub(crate) fn decompress_stream(
            &mut self,
            dst: &mut [u8],
            src: &[u8],
        ) -> Result<(usize, usize, usize), usize> {
            let mut out = ffi::ZSTD_outBuffer {
                dst: dst.as_mut_ptr().cast::<c_void>(),
                size: dst.len(),
                pos: 0,
            };
            let mut inp = ffi::ZSTD_inBuffer {
                src: src.as_ptr().cast::<c_void>(),
                size: src.len(),
                pos: 0,
            };
            let hint =
                check(unsafe { ffi::ZSTD_decompressStream(self.raw, &raw mut out, &raw mut inp) })?;
            Ok((out.pos, inp.pos, hint))
        }
    }

    impl Drop for DCtx<'_> {
        fn drop(&mut self) {
            unsafe { ffi::ZSTD_freeDCtx(self.raw) };
        }
    }

    #[must_use]
    pub(crate) fn compress_bound(src_len: usize) -> usize {
        ffi::ZSTD_compressBound(src_len)
    }

    #[must_use]
    pub(crate) fn get_error_name(code: usize) -> &'static str {
        unsafe {
            let ptr = ffi::ZSTD_getErrorName(code);
            CStr::from_ptr(ptr).to_str().unwrap_or("unknown zstd error")
        }
    }
}

/// One-shot buffer<->buffer (de)compression, matching `zstd::stream::{copy_encode,
/// copy_decode}`.
pub(crate) mod stream {
    use super::{Write, io, zstd_safe};

    pub(crate) fn copy_encode(src: &[u8], dst: &mut Vec<u8>, level: i32) -> io::Result<()> {
        let mut cctx = zstd_safe::CCtx::create();
        cctx.set_parameter(zstd_safe::CParameter::CompressionLevel(level))
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        dst.clear();
        dst.reserve(zstd_safe::compress_bound(src.len()));
        cctx.compress2(dst, src)
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        Ok(())
    }

    /// Streams one or more concatenated zstd frames from `src` into `dst`.
    ///
    /// Uses `ZSTD_decompressStream` in a loop rather than a single-shot decompress,
    /// so it works regardless of whether the frame carries an embedded content
    /// size (one-shot `compress2`-produced frames always do; frames written by a
    /// true streaming encoder often don't).
    pub(crate) fn copy_decode<W: Write>(src: &[u8], mut dst: W) -> io::Result<()> {
        const CHUNK: usize = 64 * 1024;

        let mut dctx = zstd_safe::DCtx::create();
        let mut out_buf = vec![0u8; CHUNK];
        let mut remaining = src;
        loop {
            let (written, consumed, hint) = dctx
                .decompress_stream(&mut out_buf, remaining)
                .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
            if written > 0 {
                dst.write_all(&out_buf[..written])?;
            }
            remaining = &remaining[consumed..];
            if hint == 0 {
                if remaining.is_empty() {
                    return Ok(());
                }
                // A complete frame was flushed but more input remains - could be a
                // concatenated frame, keep going.
                continue;
            }
            if consumed == 0 && written == 0 {
                return Err(io::Error::other(
                    "zstd decompression stalled (truncated input?)",
                ));
            }
        }
    }
}

/// Buffer-then-finish `Write` adapter, matching `zstd::Encoder`'s external shape.
///
/// Unlike the real crate, this is not incremental: writes accumulate into an
/// internal buffer and are only compressed and flushed to `inner` on `Drop` (or an
/// explicit [`Encoder::finish`]). Both real call sites (`vbq/index.rs`) already build
/// a complete in-memory buffer before compressing, so this is behavior-equivalent for
/// this crate's usage - it just doesn't support genuinely unbounded incremental
/// streams.
pub(crate) struct Encoder<W: Write> {
    inner: W,
    buf: Vec<u8>,
    level: i32,
    finished: bool,
}

impl<W: Write> Encoder<W> {
    /// Infallible here, but returns `io::Result` to match `zstd::Encoder::new`'s
    /// signature - both backends must expose the same call-site API.
    #[allow(clippy::unnecessary_wraps)]
    pub(crate) fn new(inner: W, level: i32) -> io::Result<Self> {
        Ok(Self {
            inner,
            buf: Vec::new(),
            level,
            finished: false,
        })
    }

    /// No-op passthrough for call-site compatibility with `zstd::Encoder::auto_finish`
    /// - this backend always finishes on `Drop`.
    #[must_use]
    pub(crate) fn auto_finish(self) -> Self {
        self
    }

    fn finish_inner(&mut self) -> io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        let mut compressed = Vec::new();
        stream::copy_encode(&self.buf, &mut compressed, self.level)?;
        self.inner.write_all(&compressed)
    }
}

impl<W: Write> Write for Encoder<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.buf.extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // Real compression + write to `inner` happens once, on `Drop` - see the
        // module-level docs on `Encoder`.
        Ok(())
    }
}

impl<W: Write> Drop for Encoder<W> {
    fn drop(&mut self) {
        // Mirrors `zstd::Encoder`'s drop-finishes behavior. Errors here can't be
        // propagated from `Drop`; this matches `AutoFinishEncoder`'s best-effort
        // semantics in the real crate.
        let _ = self.finish_inner();
    }
}

/// Eager-decode-on-first-read `Read` adapter, matching `zstd::Decoder`'s external
/// shape. See [`Encoder`] for why this is one-shot rather than incremental.
pub(crate) struct Decoder<R> {
    inner: Option<R>,
    buf: Vec<u8>,
    pos: usize,
}

impl<R: Read> Decoder<R> {
    /// Infallible here, but returns `io::Result` to match `zstd::Decoder::new`'s
    /// signature - both backends must expose the same call-site API.
    #[allow(clippy::unnecessary_wraps)]
    pub(crate) fn new(inner: R) -> io::Result<Self> {
        Ok(Self {
            inner: Some(inner),
            buf: Vec::new(),
            pos: 0,
        })
    }

    fn ensure_decoded(&mut self) -> io::Result<()> {
        if let Some(mut inner) = self.inner.take() {
            let mut compressed = Vec::new();
            inner.read_to_end(&mut compressed)?;
            stream::copy_decode(&compressed, &mut self.buf)?;
        }
        Ok(())
    }
}

impl<R: Read> Read for Decoder<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        self.ensure_decoded()?;
        let remaining = &self.buf[self.pos..];
        let n = remaining.len().min(out.len());
        out[..n].copy_from_slice(&remaining[..n]);
        self.pos += n;
        Ok(n)
    }
}
