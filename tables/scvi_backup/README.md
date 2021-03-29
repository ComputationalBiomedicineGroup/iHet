The SCVI backup is used to reproduce the cell-type clustering, as scVI is
not reproducible across different runs and different systems.

Not that for the cluster, we use the scVI coordinates stored in
the `h5ad` object. The full model is from a later run (and used
for DE analysis).
