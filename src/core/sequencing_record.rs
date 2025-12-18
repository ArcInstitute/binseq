use anyhow::{Result, bail};

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
    #[inline]
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
    #[inline]
    pub fn size(&self) -> usize {
        self.s_seq.len()
            + self.s_qual.map_or(0, |q| q.len())
            + self.s_header.map_or(0, |h| h.len())
            + self.x_seq.map_or(0, |q| q.len())
            + self.x_qual.map_or(0, |q| q.len())
            + self.x_header.map_or(0, |h| h.len())
            + self.flag.map_or(0, |f| f.to_le_bytes().len())
    }

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.x_seq.is_some()
    }
    #[inline]
    pub fn has_flags(&self) -> bool {
        self.flag.is_some()
    }
    #[inline]
    pub fn has_headers(&self) -> bool {
        self.s_header.is_some() || self.x_header.is_some()
    }
    #[inline]
    pub fn has_qualities(&self) -> bool {
        self.s_qual.is_some() || self.x_qual.is_some()
    }
}

#[derive(Default)]
pub struct SequencingRecordBuilder<'a> {
    s_seq: Option<&'a [u8]>,
    s_qual: Option<&'a [u8]>,
    s_header: Option<&'a [u8]>,
    x_seq: Option<&'a [u8]>,
    x_qual: Option<&'a [u8]>,
    x_header: Option<&'a [u8]>,
    flag: Option<u64>,
}
impl<'a> SequencingRecordBuilder<'a> {
    pub fn s_seq(mut self, s_seq: &'a [u8]) -> Self {
        self.s_seq = Some(s_seq);
        self
    }
    pub fn s_qual(mut self, s_qual: &'a [u8]) -> Self {
        self.s_qual = Some(s_qual);
        self
    }
    pub fn opt_s_qual(mut self, s_qual: Option<&'a [u8]>) -> Self {
        self.s_qual = s_qual;
        self
    }
    pub fn s_header(mut self, s_header: &'a [u8]) -> Self {
        self.s_header = Some(s_header);
        self
    }
    pub fn opt_s_header(mut self, s_header: Option<&'a [u8]>) -> Self {
        self.s_header = s_header;
        self
    }
    pub fn x_seq(mut self, x_seq: &'a [u8]) -> Self {
        self.x_seq = Some(x_seq);
        self
    }
    pub fn opt_x_seq(mut self, x_seq: Option<&'a [u8]>) -> Self {
        self.x_seq = x_seq;
        self
    }
    pub fn x_qual(mut self, x_qual: &'a [u8]) -> Self {
        self.x_qual = Some(x_qual);
        self
    }
    pub fn opt_x_qual(mut self, x_qual: Option<&'a [u8]>) -> Self {
        self.x_qual = x_qual;
        self
    }
    pub fn x_header(mut self, x_header: &'a [u8]) -> Self {
        self.x_header = Some(x_header);
        self
    }
    pub fn opt_x_header(mut self, x_header: Option<&'a [u8]>) -> Self {
        self.x_header = x_header;
        self
    }
    pub fn flag(mut self, flag: u64) -> Self {
        self.flag = Some(flag);
        self
    }
    pub fn opt_flag(mut self, flag: Option<u64>) -> Self {
        self.flag = flag;
        self
    }
}

impl<'a> SequencingRecordBuilder<'a> {
    pub fn build(self) -> Result<SequencingRecord<'a>> {
        if self.s_seq.is_none() {
            bail!("Missing s_seq on building sequencing record");
        }
        Ok(SequencingRecord {
            s_seq: self.s_seq.unwrap(),
            s_qual: self.s_qual,
            s_header: self.s_header,
            x_seq: self.x_seq,
            x_qual: self.x_qual,
            x_header: self.x_header,
            flag: self.flag,
        })
    }
}
