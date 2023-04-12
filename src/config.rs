use derive_builder::Builder;

#[derive(Debug, Clone)]
/// For the lack of a better name, a collection of user-specified "hyper-parameters" for the program
pub struct ExternalContext {
    pub hmm_size_lb: usize,
    pub show_progress: bool,
    pub io_bound: bool,
    pub trim: bool,
    pub only_queries: bool,
    pub db: Option<sled::Db>,
}
