def compute_fdr_from_labels(labels, target_fdr=5):
    """
    Compute xi-style fdr based on labels IMPORTANT: labels must be ordered by descending score.

    Parameters
    ----------
    labels : Labels in correct ordering to compute FDR

    Returns
    -------
    fdr_index: List with FDR values
    psm_hits: List with psms at a given FDR
    """
    max_fdr = 0
    decoy_decoy_count = 0
    target_decoy_count = 0
    target_target_count = 0
    fdr_index = []
    hits = []
    for l in labels:
        if l == 'DD':
            decoy_decoy_count += 1
        elif l == 'TD':
            target_decoy_count +=1
        elif l == 'TT':
            target_target_count += 1
        try:
            fdr = (target_decoy_count - decoy_decoy_count) / float(target_target_count)
            if fdr > max_fdr:
                max_fdr = fdr
            else:
                fdr = max(max_fdr, fdr)
        except:
            fdr = 0
        fdr_index.append(fdr*100)
        hits.append(target_target_count)
        if fdr*100 >= target_fdr:
            return fdr_index, hits
    return fdr_index, hits


def compute_number_of_entries_at_fdr(labels):
    """
    Compute xi-style fdr based on labels IMPORTANT: labels must be ordered by descending score.

    Parameters
    ----------
    labels : Labels in correct ordering to compute FDR

    Returns
    -------
    fdr_index: List with FDR values
    psm_hits: List with psms at a given FDR
    """
    max_fdr = 0
    decoy_decoy_count = 0
    target_decoy_count = 0
    target_target_count = 0
    fdr_index = []
    hits = []
    for l in labels:
        if l == 'DD':
            decoy_decoy_count += 1
        elif l == 'TD':
            target_decoy_count +=1
        elif l == 'TT':
            target_target_count += 1
        try:
            fdr = (target_decoy_count - decoy_decoy_count) / float(target_target_count)
            if fdr > max_fdr:
                max_fdr = fdr
            else:
                fdr = max(max_fdr, fdr)
        except:
            fdr = 0
        fdr_index.append(fdr*100)
        hits.append(target_target_count)
    return fdr_index, hits

