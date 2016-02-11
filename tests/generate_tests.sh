jf count -m 13 -s 100M -C with_ns.fasta
jf histo mer_counts.jf > with_ns.13mer_hist.canonical
jf count -m 13 -s 100M with_ns.fasta
jf histo mer_counts.jf > with_ns.13mer_hist.not_canonical
jf count -m 15 -s 100M -C with_ns.fasta
jf histo mer_counts.jf > with_ns.15mer_hist.canonical
jf count -m 15 -s 100M with_ns.fasta
jf histo mer_counts.jf > with_ns.15mer_hist.not_canonical
jf count -m 17 -s 100M -C with_ns.fasta
jf histo mer_counts.jf > with_ns.17mer_hist.canonical
jf count -m 17 -s 100M with_ns.fasta
jf histo mer_counts.jf > with_ns.17mer_hist.not_canonical

rm mer_counts.jf
