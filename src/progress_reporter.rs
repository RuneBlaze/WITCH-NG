use std::{
    sync::{
        atomic::{AtomicUsize, Ordering},
        mpsc::Receiver,
        Arc,
    },
    time::{Duration, Instant},
};
use tracing::info;

pub fn progress_reporter(
    counter: &Arc<AtomicUsize>,
    total: usize,
    interval: Duration,
    stage_name: &str,
    terminate_signal: Receiver<bool>,
) {
    let start = Instant::now();
    let mut last_report = start;
    while counter.load(Ordering::Relaxed) < total && !terminate_signal.try_recv().is_ok() {
        let now = Instant::now();
        if now.duration_since(last_report) >= interval {
            let progress = counter.load(Ordering::Relaxed);
            let percent = (progress as f64 / total as f64) * 100.0;
            info!(
                "stage: {}, work done: {:.2}% ({}/{})",
                stage_name, percent, progress, total
            );
            last_report = now;
        }
        std::thread::sleep(Duration::from_millis(300));
    }
}
