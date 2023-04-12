use crate::{
    adder::{add_queries, AdderContext},
    config::ExternalContext,
    melt::oneshot_melt,
    score_calc::ScoringCtxt,
    structures::CrucibleCtxt,
};
use anyhow::bail;
use std::{
    fs::{self, File},
    io::BufReader,
    path::PathBuf,
    time::Instant,
};
use tracing::info;

pub fn combined_analysis(
    input_path: PathBuf,
    mut backbone_path: PathBuf,
    output_path: PathBuf,
    ehmm_path: Option<PathBuf>,
    tree_path: Option<PathBuf>,
    config: &ExternalContext,
) -> anyhow::Result<()> {
    if config.trim || config.only_queries {
        bail!("Trimming and only-queries are not implemented yet");
    }
    if config.db.as_ref().map(|it| it.was_recovered()) == Some(true) && ehmm_path.is_none() {
        info!("recovered from checkpoint file, trying to reuse existing eHMM");
        backbone_path.set_extension("ehmm");
        if !backbone_path.exists() {
            bail!("checkpoint file exists but eHMM does not. Please provide the eHMM path or remove the checkpoint");
        }
    }
    // we first decide the eHMM path and also the backbone MSA path
    let (actual_backbone_path, ehmm_ctxt, ehmm_path) = if fs::metadata(&backbone_path)?.is_dir() {
        let crucible_ctxt: CrucibleCtxt =
            serde_json::from_reader(BufReader::new(File::open(backbone_path.join("melt.json"))?))?;
        let bb_path = backbone_path.join("subsets").join("0.afa");
        (bb_path, crucible_ctxt, backbone_path)
    } else {
        let actual_ehmm_dir = if let Some(ehmm_path) = ehmm_path {
            ehmm_path
        } else {
            let mut ehmm_path = backbone_path.clone();
            ehmm_path.set_extension("ehmm");
            ehmm_path
        };
        let ctxt = oneshot_melt(
            &backbone_path,
            &tree_path.expect("building eHMM must use a backbone tree"),
            &actual_ehmm_dir,
            config,
        )?;
        (backbone_path, ctxt, actual_ehmm_dir)
    };
    // then we start scoring everything
    let scorer = ScoringCtxt::from_ehmms_ctxt(ehmm_path.clone(), ehmm_ctxt, &input_path)?;
    // scoring finished
    let t = Instant::now();
    let scored = scorer.produce_payload(config)?;
    let elapsed = t.elapsed();
    info!(
        "all-against-all hmmsearch (with adjusted bitscore calculation) took {:?}",
        elapsed
    );
    let adder = AdderContext::from_scoring_ctxt(&ehmm_path, scorer, scored)?;
    add_queries(adder, &output_path, &actual_backbone_path)?;
    Ok(())
}
