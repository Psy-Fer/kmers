use std::{collections::HashMap, process::exit};
use indicatif::{ ProgressBar,  ProgressStyle};
use std::time::Instant;
use noodles::fasta;
use std::fs::File;
use std::io::BufReader;
use std::env;
use std::hash::RandomState;

#[derive(Debug, Clone)]
pub struct KmerData {
    pub count: u64,
    pub positions: Vec<u64>,
    pub scores: Vec<u64>,
    pub zeros: u64,
}

impl KmerData {
    pub fn new(pos: u64, s_len: usize) -> Self{
        let mut positions = Vec::with_capacity(s_len);
        let scores = Vec::with_capacity(s_len-1);
        positions.push(pos);
        Self { 
            count: 1,
            positions: positions,
            scores: scores,
            zeros: 0,
         }
        }
    
    pub fn update_positions(&mut self, pos: u64) {
        self.count += 1;
        self.positions.push(pos);
    }

    pub fn update_scores(&mut self, score: u64){
        self.scores.push(score);
        if score == 0 {
            self.zeros += 1;
        }
    }

    pub fn reduce(&mut self) {
        self.positions.shrink_to_fit();
        self.scores.shrink_to_fit();
    }
}




// get kmers of length k for input sequence
fn get_kmers(s: &[u8], mut h:HashMap<String, KmerData>, k: u64) -> HashMap<String, KmerData>{
    fn unsigned64tosigned64(v: u64) -> i64 {
        i64::try_from(v).ok().unwrap()
    }
    let l: u64 = s.len().try_into().unwrap(); // sequence length
    let r: u64 = (l-k+1).try_into().unwrap(); // range length of sequence len - kmer len +1
    // loop over the sequence 1 step at a time and take a slice the length of k and store in HashMap
    let style2 = ProgressStyle::with_template("[{elapsed_precise}] {bar:80.orange/yellow} {pos:>7}/{len:7} kmer_processing {per_sec} {eta_precise} {msg}")
        .unwrap();
    let pb2: ProgressBar = ProgressBar::new((0..r as usize).len() as u64).with_style(style2);
    for i in pb2.wrap_iter(0..r) {
        let vkmer: Vec<u8> = s[i as usize..i as usize+k as usize].try_into().unwrap();
        let skmer: String = String::from_utf8(vkmer).expect("Found invalid UTF-8");
        // check if string contains N, if so, skip
        if skmer.contains("N") {
            continue;
        }
        // if i >= 15_000_000 {
        //     exit(0)
        // }
        // if kmer exists, update the count and position information, else insert a new entry
        // I could just intervene when i get an overlap position, and do the whole logic i do later to eliminate negative scores...
        if h.contains_key(&skmer) {
            // get a mutable reference to the struct, so changes don't have to be pushed back in
            let data: &mut KmerData = h.get_mut(&skmer).expect("kmer not found, but it should have been...?");
            // do score check...if negative, it's an overlapping kmer, ignore.
            // get the last element in position, that is the prev value (we switching kmers a lot here, so need to do this)
            let prev_score: u64 = data.positions.last().unwrap().clone();
            // type casting u64 into i64, because kmer value is never gonna be that big
            let po: i64 = unsigned64tosigned64(i);
            let pr: i64 = unsigned64tosigned64(prev_score);
            let k_len: i64 = unsigned64tosigned64(k);
            let score: i64 = po - pr - k_len;
            if score < 0{
                // skip this kmer
                continue;
            }
            else{
                //add this kmer position and score
                data.update_scores(score as u64);
                data.update_positions(i);
                // KmerData::update_scores(data, score as u64);
                // KmerData::update_positions(data, i);
            }

            // remove this clone
            // let new_data = data.clone();
            // h.insert(skmer, new_data);
        }
        else {
            let data: KmerData = KmerData::new(i, s.len());
            h.insert(skmer, data);
        }
    }
    pb2.abandon_with_message("✅");
    h.shrink_to_fit();
    h
}

// get a range of kmers up to max_kmer length for input sequence
fn get_kmer_range(min_kmer: u64, s: &[u8], mut kc: HashMap<u64, HashMap<String,KmerData>>, r: u64) -> HashMap<u64, HashMap<String,KmerData>>{
    let mut mrange: u64 = r + 1;
    // check the sequence is longer than the max_kmer length
    if s.len() <= r as usize{
        mrange = s.len() as u64 + 1;
    }
    let style1 = ProgressStyle::with_template("[{elapsed_precise}] {bar:80.red/purple} {pos:>7}/{len:7} get_kmers {msg}")
        .unwrap();
    let pb1 = ProgressBar::new((min_kmer as u32..mrange as u32).len() as u64).with_style(style1);
    // This should be threaded
    for k in pb1.wrap_iter(min_kmer as u32..mrange as u32)  {
        let hasher = RandomState::new();
        // max capacity is 4^k
        let kmers: HashMap<String, KmerData> = HashMap::with_capacity_and_hasher(4_usize.pow(k), hasher);
        kc.insert(k as u64, get_kmers(s, kmers, k as u64));
        // if k == mrange as u32 {
        // }
    }
    pb1.abandon_with_message("✅");
    
    kc
}



fn main() {

    let args: Vec<String> = env::args().collect();
    println!("args: {:?}", args);
    let filename = args[1].clone();
    // Read fasta file
    let mut reader = File::open(filename)
                                                .map(BufReader::new)
                                                .map(fasta::Reader::new).expect("something fucked up");

    // for each record in the fasta file, analyse and write results
    for result in reader.records() {
        let record = result.expect("invalid record");
        // first white space sep text in header
        let name = String::from_utf8(record.name().to_vec()).expect("Found invalid UTF-8");
        let seq = record.sequence();
        let seq_length = seq.len();
        // eprintln!("record: {:?}", record);
        // eprintln!("header: {}", description);
        eprintln!("name: {}", name);
        eprintln!("seq_len: {}", seq_length);
        let raw_seq: &[u8] = &seq[..];
        eprintln!("converting to upper case");
        // convert sequence to upper case
        let seq_string = String::from_utf8(raw_seq.to_vec()).expect("Found invalid UTF-8").to_uppercase();
        // eprintln!("seq-upper: {:?}", seq_string);
        
        eprintln!("converting to iseq");
        // convert upper-string to u8
        let i_seq: &[u8] = seq_string.as_bytes();
        // let i_seq: &[u8] = &seq[..];
        // let i_qual: Option<&[u8]>;

        let min_kmer: u64 = 2;
        let max_kmer: u64 = 6;
        
        eprintln!("Getting kmer-range...");
        let hasher = RandomState::new();
        let mut kmer_collection: HashMap<u64,HashMap<String, KmerData>> = HashMap::with_capacity_and_hasher((max_kmer-min_kmer) as usize, hasher);

        let now = Instant::now();
        kmer_collection = get_kmer_range(min_kmer, i_seq, kmer_collection, max_kmer);
        eprintln!("get_kmer_range time: {:.2?}", now.elapsed());

        // for (i, d) in kmer_collection {
        //     println!("kmer length: {}", i);
        //     for (k, j) in d {
        //         println!("kmer: {}", k);
        //         println!("KmerData: {:?}", j);
        //     }
        // }
    }

}
