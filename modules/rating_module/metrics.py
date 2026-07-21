import statistics

def calculate_n50(lengths):
    """Calculates the N50 metric for a list of sequence lengths."""
    if not lengths:
        return 0
    lengths.sort(reverse=True)
    total_length = sum(lengths)
    target = total_length / 2.0
    cum_sum = 0
    for length in lengths:
        cum_sum += length
        if cum_sum >= target:
            return length
    return 0

def calculate_basic_metrics(seq_info):
    """
    Calculates Yield Profile, Class-Specific Lengths, and Payload Score.
    Automatically detects TE classes without hardcoding.
    """
    total_seqs = len(seq_info)
    total_bases = 0
    known_te_bases = 0
    
    class_counts = {}
    length_pools = {}
    
    for seq_id, info in seq_info.items():
        norm_class = info['norm_class']
        length = info['length']
        
        total_bases += length
        
        # Calculate bases for known TEs (excluding Unknown)
        if norm_class.lower() != "unknown":
            known_te_bases += length
            
        # Dynamically build category dictionaries
        if norm_class not in class_counts:
            class_counts[norm_class] = 0
            length_pools[norm_class] = []
            
        class_counts[norm_class] += 1
        length_pools[norm_class].append(length)

    results = {
        "Yield": {
            "Total_Sequences": total_seqs,
            "Total_Bases": total_bases,
            "Distribution": {}
        },
        "Lengths": {},
        "Payload_Score": 0.0
    }

    # Calculate Payload Score (percentage of total bases belonging to known TEs)
    if total_bases > 0:
        results["Payload_Score"] = round((known_te_bases / total_bases) * 100, 2)

    # Calculate percentages for Yield Profile
    for cls, count in class_counts.items():
        pct = (count / total_seqs * 100) if total_seqs > 0 else 0
        results["Yield"]["Distribution"][cls] = {"count": count, "percent": round(pct, 2)}

    # Calculate N50 and Median for Length Metrics
    for cls, lengths in length_pools.items():
        median_val = statistics.median(lengths) if lengths else 0
        n50_val = calculate_n50(lengths)
        results["Lengths"][cls] = {
            "Median": round(median_val, 2),
            "N50": n50_val,
            "Pool_Size": len(lengths)
        }

    return results
