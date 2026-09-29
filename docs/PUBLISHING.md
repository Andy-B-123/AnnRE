# Publish the AnnRE MVP

The target repository is **https://github.com/Andy-B-123/AnnRE**. The prepared folder is a standalone source tree; it has not been published by this preparation step.

## 1. Review the release

Read README.md, ACKNOWLEDGEMENTS.md and docs/BSF_RESULTS.md. Confirm the contributor list in CITATION.cff and your institution's requirements. Choose the software license with the appropriate project owners and add its complete text as LICENSE; then update LICENSE_STATUS.md. Until then this is a research preview with no open-source license granted.

Keep raw sequencing data, sample-identifying metadata, email, slides, credentials and private project files outside this repository. The supplied examples use generic paths and synthetic sample names. Inspect your final diff before publishing.

## 2. Simplest route: GitHub Desktop

1. Download the prepared source ZIP to your computer and extract it.
2. In GitHub Desktop, sign in as Andy-B-123. Select **File → New repository**, name it **AnnRE**, and choose a local parent folder. Avoid selecting an automatic license until you have chosen one.
3. Copy the **contents** of the extracted AnnRE folder into the new repository folder. README.md and pyproject.toml should be at the repository root, not inside a second AnnRE directory. Include the hidden `.github` and `.gitignore` files.
4. Review the Changes tab. Commit with a message such as `Initial AnnRE research alpha`.
5. Click **Publish repository**. Confirm the name and account. To make the conference link accessible to everyone, clear **Keep this code private** before publishing.
6. Open https://github.com/Andy-B-123/AnnRE in a signed-out browser. Check the README graphics and Actions test results.

If the repository already exists, clone it with GitHub Desktop and copy the prepared files into that clone, then commit and push. Do not overwrite existing work without reviewing the differences.

Official guide: [Creating your first repository with GitHub Desktop](https://docs.github.com/en/desktop/overview/creating-your-first-repository-using-github-desktop).

## 3. Command-line alternative

Create an **empty** AnnRE repository under Andy-B-123 on GitHub, with no generated README, license or gitignore. In a shell, enter the prepared standalone AnnRE folder, then run:

```bash
git init -b main
git add .
git status
git diff --cached --stat
git commit -m "Initial AnnRE research alpha"
git remote add origin https://github.com/Andy-B-123/AnnRE.git
git push -u origin main
```

Authenticate using your configured GitHub credential manager or SSH setup; never put a token into these files or the remote URL. If the remote has existing commits, clone it first and copy the source tree into that clone instead of using this empty-repository procedure. Do not force-push over existing history.

Official guide: [Adding locally hosted code to GitHub](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github).

## 4. Optional splash website

The repository already includes a static website in `docs/index.html`.

On GitHub, open **Settings → Pages**. Set the source to **Deploy from a branch**, select **main**, choose **/docs**, and save. After deployment, the expected URL is https://andy-b-123.github.io/AnnRE/ . Check the actual URL reported by GitHub before sharing it. The conference QR supplied here points to the repository, so Pages is optional.

Official guide: [Configuring a GitHub Pages publishing source](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

## 5. Mark the release as an alpha

Once the tests pass and metadata/license are ready, use **Releases → Draft a new release**. Create tag `v0.1.0a1`, title `AnnRE 0.1 alpha — evidence-to-review MVP`, and select **Set as a pre-release**. Describe the available modules, ten-gene portable demonstration, separate genome-wide BSF census and limitations. No PyPI upload is needed for this MVP.

For subsequent changes, edit the local source, run the tests, review the diff, commit and push. Treat biological adjudication and independent benchmarking as release work, not as conclusions implied by the alpha label.
