use std::io;

use anyhow::Result;
use zstd::zstd_safe;

pub(crate) fn sized_compress(
    dst: &mut Vec<u8>,
    src: &[u8],
    level: u64,
    cctx: &mut zstd_safe::CCtx,
) -> Result<()> {
    // determine the maximum compressed size
    let max_z_size = zstd_safe::compress_bound(src.len());

    // resize the destination vector to the maximum compressed size
    //
    // Note: this uses uninitialized memory, but is safe because we immediately
    // follow it with a call to `compress` which overwrites the buffer.
    resize_uninit(dst, max_z_size);

    // Compress the data using the provided compression context
    let true_size = cctx
        .compress(dst, src, level as i32)
        .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;

    // resize to the true size - clipping all remaining uninitialized memory
    dst.truncate(true_size);

    Ok(())
}

pub(crate) fn extension_read<R: io::Read>(
    reader: &mut R,
    dst: &mut Vec<u8>,
    size: usize,
) -> Result<()> {
    dst.resize(size, 0);
    reader.read_exact(dst)?;
    Ok(())
}

pub(crate) fn slice_and_increment<'a>(offset: &mut usize, len: u64, bytes: &'a [u8]) -> &'a [u8] {
    let slice = &bytes[*offset..*offset + len as usize];
    *offset += len as usize;
    slice
}

/// Resize a vector to the target length without initializing new elements.
///
/// # Safety
/// The caller must ensure that all elements in the range [old_len..new_len]
/// are initialized before reading them. This is safe when immediately followed
/// by operations that write to the entire buffer (e.g., decompression).
#[inline]
pub(crate) fn resize_uninit<T>(vec: &mut Vec<T>, new_len: usize) {
    match new_len.cmp(&vec.len()) {
        std::cmp::Ordering::Greater => {
            // Growing: reserve and set length (unsafe but fast)
            vec.reserve(new_len - vec.len());
            unsafe {
                vec.set_len(new_len);
            }
        }
        std::cmp::Ordering::Less => {
            // Shrinking: truncate (safe and fast)
            vec.truncate(new_len);
        }
        std::cmp::Ordering::Equal => {
            // Same size: do nothing
        }
    }
}

pub(crate) fn calculate_offsets(values: &[u64], offsets: &mut Vec<u64>) {
    offsets.clear();
    offsets.push(0);
    for i in 1..values.len() {
        offsets.push(offsets[i - 1] + values[i - 1]);
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Span {
    offset: usize,
    length: usize,
}
impl Span {
    pub fn new(offset: usize, length: usize) -> Self {
        Span { offset, length }
    }

    pub fn new_u64(offset: u64, length: u64) -> Self {
        Span::new(offset as usize, length as usize)
    }

    pub fn range(&self) -> std::ops::Range<usize> {
        self.offset..self.offset + self.length
    }

    pub fn len(&self) -> usize {
        self.length
    }
}
