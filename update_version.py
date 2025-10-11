# update_citation.py
import tomli
import tomli_w
from oscilate.__version__ import __version__

# Update citation.txt
citation_content = f"""Vincent Mahé. (2025). OSCILATE (Version {__version__}) [Computer software]. GitHub. https://github.com/VinceECN/OSCILATE"""
with open("citation.txt", "w") as f:
    f.write(citation_content)

# Update CITATION.bib
bib_content = f"""@software{{Mahe_OSCILATE_2025,
  author = {{Mahé, Vincent}},
  title = {{OSCILATE: A Python implementation of the Method of Multiple Scales using SymPy for symbolic computation}},
  year = {{2025}},
  publisher = {{GitHub}},
  url = {{https://github.com/VinceECN/OSCILATE}},
  version = {{{__version__}}}
}}
"""
with open("CITATION.bib", "w") as f:
    f.write(bib_content)

# Update pyproject.toml
import tomli_w
with open("pyproject.toml", "rb") as f:
    toml_dict = tomli.load(f)
toml_dict["project"]["version"] = __version__
with open("pyproject.toml", "wb") as f:
    tomli_w.dump(toml_dict, f)