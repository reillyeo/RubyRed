MIMt_16s_311025_species classifier retrained on the 341F-1390R amplicon region.

Built 2026-09-30 (QIIME 2 2024.5, scikit-learn 1.4.2):
  qiime feature-classifier extract-reads \
    --i-sequences ../MIMt_16s_311025_species/MIMt_16s_311025_species_seqs.qza \
    --p-f-primer CCTACGGGNGGCWGCAG --p-r-primer ACGGGCGGTGTGTRCA \
    --p-min-length 800 --p-max-length 1300 --p-read-orientation both \
    --o-reads MIMt_16s_311025_species_amp341F1390R_seqs.qza
  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads MIMt_16s_311025_species_amp341F1390R_seqs.qza \
    --i-reference-taxonomy ../MIMt_16s_311025_species/MIMt_16s_311025_species_tax.qza \
    --o-classifier MIMt_16s_311025_species_amp341F1390R_classifier.qza

30,239 of 31,436 references contained both primer sites (amplicon 984-1040 bp);
668 of 24,441 species lost all references (listed in lost_species.txt).
Taxonomy is the species-level MIMt_16s_311025_species taxonomy (symlinked).

Intended for reads trimmed to the same primers (RubyRed_trim). Use with:
  RubyRed_trim -c <this dir>/MIMt_16s_311025_species_amp341F1390R_classifier.qza
