"""Helper functions to fetch data from NSRR datasets."""
import glob


def get_sub_visit_hyp(eeg_file, dataset, hypno_dir):
    """Extract the subject ID, visit ID and hypno_file."""
    if dataset == "mesa":
        sub = eeg_file.split("-")[2][:-4]
        visit = "visit1"
        hypno_file = hypno_dir + 'mesa-sleep-' + str(sub) + '-profusion.xml'
    elif dataset == "cfs":
        sub = eeg_file.split("-")[2][:-4]
        visit = "visit1"
        hypno_file = hypno_dir + 'cfs-visit5-' + str(sub) + '-profusion.xml'
    return sub, visit, hypno_file


def get_all_edfs(dataset, root_dir):
    """List all EDFs files for a given study."""
    if dataset == "mesa":
        eeg_dir = root_dir + 'mesa/polysomnography/edfs/'
        hypno_dir = root_dir + 'mesa/polysomnography/annotations-events-profusion/'
        all_edfs = sorted(glob.glob(eeg_dir + "*.edf"))
    elif dataset == "cfs":
        eeg_dir = root_dir + 'cfs/polysomnography/edfs/'
        hypno_dir = root_dir + 'cfs/polysomnography/annotations-events-profusion/'
        all_edfs = sorted(glob.glob(eeg_dir + "*.edf"))
    return all_edfs, hypno_dir
