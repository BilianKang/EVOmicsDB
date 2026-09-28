# Third-party material

MIT covers the authors' software and documentation only. Reused experimental data, derived source data, miRNA/lncRNA annotations, KEGG and MSigDB resources remain subject to their original access and reuse conditions. No new data licence is assigned by this packaging step.

KEGG REST captures and four KEGG-derived reference tables are excluded from the public GitHub code package and the public Zenodo analysis archive. Their original source URLs, capture identifiers and checksums remain in the snapshot manifest. KEGG-synonym descriptions have also been removed from the retained author metabolite mapping. That mapping retains the original input names, KEGG IDs and author mapping rules; its audit retains mapping decisions. This package does not grant a right to download or redistribute upstream KEGG data. Use the applicable upstream terms and any required permission: https://www.kegg.jp/kegg/legal.html and https://www.kegg.jp/kegg/rest/. See docs/REPRODUCING.md#external-annotations.

The retained gene-symbol mapping derives from org.Hs.eg.db 3.22.0, with its source/version recorded in the snapshot manifest. Its upstream package licence and database provenance apply; it is not relabelled as original MIT-licensed author data. See https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html.

Reference target databases and full experimental matrices remain external inputs; data/input_assets.json records their identities. Historical validation documents are retained in the paired full analysis archive. They do not imply that excluded reference tables are distributed in this code ZIP.
