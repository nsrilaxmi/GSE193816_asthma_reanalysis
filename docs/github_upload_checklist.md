# GitHub Upload Checklist

Before pushing:

```bash
git status --short
git add .github CITATION.cff README.md data docs environment.yml notebooks reports requirements.txt results scripts workflow
git commit -m "Add Snakemake workflow for GSE193816 asthma reanalysis"
```

Create a new empty GitHub repository, then:

```bash
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/GSE193816_asthma_reanalysis.git
git push -u origin main
```

After pushing:

- Confirm GitHub Actions completes successfully.
- Confirm README figures render correctly on GitHub.
- Confirm the optional Snakemake instructions in `workflow/README.md` are visible.
- Keep large `.h5ad` files and generated result tables out of Git unless there is a specific release reason.
