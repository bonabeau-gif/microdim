# Publishing on GitHub

The root of this directory is the repository root. The repository contains no development diary, change log, or unpublished experiment chronology.

A new public GitHub repository can be created and populated with:

```bash
git init -b main
git add .
git commit -m "Add manuscript reproducibility repository"
gh repo create gut-resource-dimension --public --source=. --remote=origin --push
```

The commands require Git, the GitHub CLI, and prior `gh auth login` on the user's machine. Repository creation and the first upload must be authorized in the owner's GitHub account.

The `.github/workflows/checks.yml` workflow verifies unit tests, reference results, and figure generation on each push. For an archival DOI, enable GitHub-to-Zenodo archival after repository publication and place the DOI in the manuscript's Code Availability statement.

Third-party model files downloaded into `data/third_party/` are ignored by Git. Cite their original publications and obey their respective license terms.
