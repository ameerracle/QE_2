from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).resolve().parent


setup(
    name="qe2-reader",
    version="0.1.0",
    description="Utilities for reading Quantum ESPRESSO output files",
    long_description=(ROOT / "README.md").read_text() if (ROOT / "README.md").exists() else "",
    long_description_content_type="text/markdown",
    py_modules=["readqe"],
    packages=find_packages(),
    install_requires=["ase", "numpy"],
    entry_points={
        "console_scripts": [
            "readqe=readqe:_main",
        ]
    },
    python_requires=">=3.10",
)