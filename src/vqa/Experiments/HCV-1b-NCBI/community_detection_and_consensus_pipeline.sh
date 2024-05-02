


python ../../make_graph.py --results predictions.tsv --output communities.tsv

python ../../make_consensus.py --communities communities.tsv --output consensus.tsv

python ../../evaluate_consensus.py --consensus consensus.tsv --output consensus_evaluation.tsv --gt_dir 

