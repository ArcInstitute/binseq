#[derive(Clone, Copy, Default)]
pub struct SequencingRecord<'a> {
    pub(crate) s_seq: &'a [u8],
    pub(crate) s_qual: Option<&'a [u8]>,
    pub(crate) s_header: Option<&'a [u8]>,
    pub(crate) x_seq: Option<&'a [u8]>,
    pub(crate) x_qual: Option<&'a [u8]>,
    pub(crate) x_header: Option<&'a [u8]>,
    pub(crate) flag: Option<u64>,
}
impl<'a> SequencingRecord<'a> {
    pub fn new(
        s_seq: &'a [u8],
        s_qual: Option<&'a [u8]>,
        s_header: Option<&'a [u8]>,
        x_seq: Option<&'a [u8]>,
        x_qual: Option<&'a [u8]>,
        x_header: Option<&'a [u8]>,
        flag: Option<u64>,
    ) -> Self {
        Self {
            s_seq,
            s_qual,
            s_header,
            x_seq,
            x_qual,
            x_header,
            flag,
        }
    }

    /// Returns the size of the record in bytes
    pub fn size(&self) -> usize {
        self.s_seq.len()
            + self.s_qual.map_or(0, |q| q.len())
            + self.s_header.map_or(0, |h| h.len())
            + self.x_seq.map_or(0, |q| q.len())
            + self.x_qual.map_or(0, |q| q.len())
            + self.x_header.map_or(0, |h| h.len())
            + self.flag.map_or(0, |f| f.to_le_bytes().len())
    }

    pub fn is_paired(&self) -> bool {
        self.x_seq.is_some()
    }
}
