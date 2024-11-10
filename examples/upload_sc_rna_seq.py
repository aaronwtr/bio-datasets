"""We upload asymmetric units.

Ultimately what we want to be able to do is to infer the assembly from the coordinates for a single repeating unit.

Before running this script, download the single-cell data.

Download from https://zenodo.org/records/7041849
"""
import argparse
import glob
import os
import shutil
import subprocess
import tempfile
import scanpy as sc

from bio_datasets import Dataset, Features, NamedSplit, Value
from bio_datasets.features import AtomArrayFeature, StructureFeature


def examples_generator(
        sc_data_path
):
    adata = sc.read_h5ad(sc_data_path)
    for gene_name in adata.var_names:
        yield {
            "id": gene_name,
            "expression": adata[:, gene_name].X,
        }


def main(args):
    features = Features(
        id=Value("string"),
        structure=AtomArrayFeature() if args.as_array else StructureFeature(),
    )

    with tempfile.TemporaryDirectory(dir=args.temp_dir) as temp_dir:
        ds = Dataset.from_generator(
            examples_generator,
            gen_kwargs={
                "pair_codes": args.pair_codes,
                "pdb_download_dir": args.pdb_download_dir,
                "compress": args.compress,
                "remove_cif": args.remove_cif,
            },
            features=features,
            cache_dir=temp_dir,
            split=NamedSplit("train"),
            num_proc=args.num_proc,
        )
        ds.push_to_hub(
            "biodatasets/pdb",
            config_name=args.config_name or "default",
            max_shard_size="350MB",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_name", type=str, default=None)
    parser.add_argument(
        "--pair_codes", nargs="+", help="PDB 2-letter codes", default=None
    )
    parser.add_argument(
        "--backbone_only", action="store_true", help="Whether to drop sidechains"
    )
    parser.add_argument(
        "--as_array", action="store_true", help="Whether to return an array"
    )
    parser.add_argument(
        "--pdb_download_dir",
        type=str,
        default="data/pdb",
        help="Directory to download PDBs to",
    )
    parser.add_argument(
        "--temp_dir",
        type=str,
        default=None,
        help="Temporary directory (for caching built dataset)",
    )
    parser.add_argument(
        "--compress",
        action="store_true",
        help="Whether to compress the compressed bcif with gzip",
    )
    parser.add_argument(
        "--remove_cif",
        action="store_true",
        help="Whether to remove the original CIF files after conversion",
    )
    parser.add_argument(
        "--num_proc",
        type=int,
        default=1,
        help="Number of processes to use",
    )
    args = parser.parse_args()
    os.makedirs(args.pdb_download_dir, exist_ok=True)
    if args.temp_dir is not None:
        os.makedirs(args.temp_dir, exist_ok=True)

    main(args)
