# GitHub Upload Checklist

Before pushing:

```bash
git status --short
git add .
git commit -m "Initial GSE193816 asthma reanalysis workflow"
```

Create a new empty GitHub repository, then:

```bash
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/GSE193816_asthma_reanalysis.git
git push -u origin main
```

After pushing:

- Replace `YOUR_USERNAME` in `CITATION.cff`.
- Replace `Add your name` in `CITATION.cff`.
- Add selected final figures to `docs/figures/` only after reviewing them.
- Keep large `.h5ad` files and generated result tables out of Git unless there is a specific release reason.

