# Alpha preparation checks

The AnnRE 0.1.0a1 preparation passed:

- Nine synthetic Python unit/integration tests, including an end-to-end prepare/extract/report workflow.
- Wheel build and installation into an isolated temporary target; installed module reports version `0.1.0a1`.
- Local Markdown file-link checks.
- Independent QR decoding of both the standalone QR PNG and composed conference slide to `https://github.com/Andy-B-123/AnnRE`.
- Visual inspection of the local Firefox splash screenshot and 16:9 conference graphic.

The AnnRE source was adapted from the existing tested prototype. This branding/release preparation did not rerun the BSF dataset. The earlier ten-gene demonstration and project-specific genome-wide census are described separately. GitHub Actions is configured for Python 3.10 and 3.12; those remote CI jobs will run after the repository is uploaded.

These checks verify implementation and presentation behavior, not biological accuracy. This preparation does not publish to GitHub, select an open-source license or establish a complete scientific author list.
