use rayon::{ThreadPool, ThreadPoolBuilder};

#[derive(Debug, Clone)]
/// For the lack of a better name, a collection of user-specified "hyper-parameters" for the program
pub struct ExternalContext {
    pub hmm_size_lb: usize,
    pub show_progress: bool,
    pub io_bound: bool,
    pub trim: bool,
    pub only_queries: bool,
    pub num_workers: usize,
    pub num_threads_per_worker: usize,
    pub db: Option<sled::Db>,
}

impl ExternalContext {
    pub fn total_threads(&self) -> usize {
        self.num_workers * self.num_threads_per_worker
    }

    pub fn create_full_pool(&self) -> ThreadPool {
        ThreadPoolBuilder::new()
            .num_threads(self.total_threads())
            .build()
            .expect("Failed to create thread pool.")
    }
}

pub fn create_pool(num_threads: usize) -> ThreadPool {
    ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to create thread pool.")
}
