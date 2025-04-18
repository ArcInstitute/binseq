use std::arch::aarch64::*;

/// **AArch64 NEON availability** – always `true` on ARM64.
#[inline]
pub fn is_neon_supported() -> bool {
    true
}

/// Encode 16 ASCII nucleotides (`A`, `C`, `G`, `T`) into a single `u32`.
///
/// Output layout: nt0 → bits 0‑1 … nt15 → bits 30‑31 (little‑endian).
#[inline(always)]
pub unsafe fn encode_16_nucleotides(nucs: uint8x16_t) -> u32 {
    // 1. ASCII → 2‑bit codes: code = ((b >> 1) ^ (b >> 2)) & 3
    let t1 = vshrq_n_u8(nucs, 1);
    let t2 = vshrq_n_u8(nucs, 2);
    let code = vandq_u8(veorq_u8(t1, t2), vdupq_n_u8(3));

    // 2. Pack two codes into one 4‑bit nibble
    let even = vuzp1q_u8(code, code); // c0, c2, …, c14
    let odd = vuzp2q_u8(code, code); // c1, c3, …, c15
    let nibbles = vorrq_u8(even, vshlq_n_u8(odd, 2));

    // 3. Pack two nibbles into one byte
    let even_b = vuzp1q_u8(nibbles, nibbles); // p0, p2, p4, p6
    let odd_b = vuzp2q_u8(nibbles, nibbles); // p1, p3, p5, p7
    let packed = vorrq_u8(even_b, vshlq_n_u8(odd_b, 4));

    // 4. Return the first lane (lower 32 bits)
    vgetq_lane_u32(vreinterpretq_u32_u8(packed), 0)
}

/// Return `true` if every byte in `v` is a valid nucleotide (case‑insensitive).
#[inline(always)]
unsafe fn valid_block(v: uint8x16_t) -> bool {
    let lower = vorrq_u8(v, vdupq_n_u8(0x20));
    let is_a = vceqq_u8(lower, vdupq_n_u8(b'a'));
    let is_c = vceqq_u8(lower, vdupq_n_u8(b'c'));
    let is_g = vceqq_u8(lower, vdupq_n_u8(b'g'));
    let is_t = vceqq_u8(lower, vdupq_n_u8(b't'));
    let ok = vorrq_u8(is_a, vorrq_u8(is_c, vorrq_u8(is_g, is_t)));
    vminvq_u8(ok) == 0xFF
}

/// Encode an arbitrary‑length ASCII slice into packed 2‑bit words (`u64`).
///
/// * 32 nt per word.
/// * `output` must be large enough; otherwise `Err(())` is returned.
/// * On any invalid byte the function zero‑fills `output` and returns `Err(())`.
#[cfg(target_arch = "aarch64")]
pub unsafe fn encode_nucleotides_simd(input: &[u8], output: &mut [u64]) -> Result<(), ()> {
    let need = (input.len() + 31) / 32;
    if output.len() < need {
        return Err(());
    }

    output.fill(0);

    let mut ip = input.as_ptr();
    let mut left = input.len();
    let mut out = output.as_mut_ptr();

    // Vector loop: 32 nt → 1 u64
    while left >= 32 {
        let v0 = vld1q_u8(ip);
        let v1 = vld1q_u8(ip.add(16));
        if !valid_block(v0) || !valid_block(v1) {
            return Err(());
        }
        *out = (encode_16_nucleotides(v0) as u64) | ((encode_16_nucleotides(v1) as u64) << 32);

        ip = ip.add(32);
        left -= 32;
        out = out.add(1);
    }

    // Scalar tail (≤ 31 nt)
    if left != 0 {
        let mut tail = 0u64;
        for i in 0..left {
            tail |= match *ip.add(i) | 0x20 {
                b'a' => 0u64,
                b'c' => 1u64,
                b'g' => 2u64,
                b't' => 3u64,
                _ => return Err(()),
            } << (2 * i);
        }
        *out = tail;
    }
    Ok(())
}

/// Decode 16 packed 2‑bit codes (`u32`) to ASCII (`A`, `C`, `G`, `T`).
#[inline(always)]
pub unsafe fn decode_16_nucleotides(encoded: u32, dst: *mut u8) {
    // 1. Broadcast the word to four lanes
    let val = vdupq_n_u32(encoded);
    let mask = vdupq_n_u32(3);

    // 2. Extract 2‑bit fields (negative counts = right shift)
    #[inline(always)]
    const fn shv(a: i32, b: i32, c: i32, d: i32) -> int32x4_t {
        unsafe { core::mem::transmute([a, b, c, d]) }
    }
    let c0 = vandq_u32(vshlq_u32(val, shv(0, -2, -4, -6)), mask);
    let c1 = vandq_u32(vshlq_u32(val, shv(-8, -10, -12, -14)), mask);
    let c2 = vandq_u32(vshlq_u32(val, shv(-16, -18, -20, -22)), mask);
    let c3 = vandq_u32(vshlq_u32(val, shv(-24, -26, -28, -30)), mask);

    // 3. Narrow u32 → u8 and assemble 16 indices
    let idx: uint8x16_t = vcombine_u8(
        vmovn_u16(vcombine_u16(vmovn_u32(c0), vmovn_u32(c1))),
        vmovn_u16(vcombine_u16(vmovn_u32(c2), vmovn_u32(c3))),
    );

    // 4. LUT: 0→A, 1→C, 2→G, 3→T
    let lut: uint8x16_t = core::mem::transmute([
        b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T', b'A', b'C', b'G',
        b'T',
    ]);
    let ascii = vqtbl1q_u8(lut, idx);

    // 5. Store
    vst1q_u8(dst, ascii);
}

/// Decode a packed 2‑bit stream (`u64` words) back to ASCII nucleotides.
#[cfg(target_arch = "aarch64")]
pub unsafe fn decode_nucleotides_simd(
    input: &[u64],
    len: usize,
    output: &mut [u8],
) -> Result<(), ()> {
    if len > output.len() {
        return Err(());
    }

    let chunk = 32;
    let chunks = len / chunk;

    for i in 0..chunks {
        let w = input.get(i).copied().unwrap_or(0);
        decode_16_nucleotides(w as u32, output.as_mut_ptr().add(i * chunk));
        decode_16_nucleotides((w >> 32) as u32, output.as_mut_ptr().add(i * chunk + 16));
    }

    // Scalar tail
    let lut = [b'A', b'C', b'G', b'T'];
    for j in (chunks * chunk)..len {
        let idx = ((input[j / 32] >> (2 * (j % 32))) & 3) as usize;
        output[j] = lut[idx];
    }
    Ok(())
}

// Convenience wrappers ---------------------------------------------------

pub fn encode(seq: &[u8], out: &mut Vec<u64>) -> Result<(), ()> {
    out.resize((seq.len() + 31) / 32, 0);
    unsafe { encode_nucleotides_simd(seq, out) }
}

pub fn decode(enc: &[u64], len: usize, out: &mut Vec<u8>) -> Result<(), ()> {
    out.resize(len, 0);
    unsafe { decode_nucleotides_simd(enc, len, out) }
}

// Tests ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neon_available() {
        assert!(is_neon_supported());
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn roundtrip_16() {
        let src = b"ACGTACGTACGTACGT";
        let mut enc = [0u64; 1];
        let mut dec = [0u8; 16];
        unsafe {
            encode_nucleotides_simd(src, &mut enc).unwrap();
            decode_nucleotides_simd(&enc, src.len(), &mut dec).unwrap();
        }
        assert_eq!(src, &dec);
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn reject_invalid() {
        let src = b"ACGTACGTNACGTACGT"; // contains 'N'
        let mut enc = [0u64; 1];
        assert!(unsafe { encode_nucleotides_simd(src, &mut enc) }.is_err());
    }
}
