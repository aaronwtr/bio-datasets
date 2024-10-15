from typing import List

import bio_datasets
import datasets

from datasets.packaged_modules.folder_based_builder import folder_based_builder


logger = datasets.utils.logging.get_logger(__name__)


class StructureFolderConfig(folder_based_builder.FolderBasedBuilderConfig):
    """BuilderConfig for StructureFolder."""

    drop_labels: bool = None
    drop_metadata: bool = None

    def __post_init__(self):
        super().__post_init__()


class ProteinStructureFolder(folder_based_builder.FolderBasedBuilder):
    BASE_FEATURE = bio_datasets.ProteinStructureFeature
    BASE_COLUMN_NAME = "structure"
    BUILDER_CONFIG_CLASS = StructureFolderConfig
    EXTENSIONS: List[str]  # definition at the bottom of the script


STRUCTURE_EXTENSIONS = [
    ".fcz",
    ".pdb",
    ".cif",
]
ProteinStructureFolder.EXTENSIONS = STRUCTURE_EXTENSIONS
