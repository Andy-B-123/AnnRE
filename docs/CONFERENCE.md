# Conference assets

- `assets/annre-talk-slide.png` / `.svg` / `.pdf`: 16:9 conference slide with QR and attribution.
- `assets/annre-banner.svg` / `.png`: README branding.
- `assets/annre-overview.svg` / `.png`: editable vector and slide-ready workflow.
- `assets/annre-splash.png`: screenshot of the local HTML splash page, not a screenshot of a live GitHub repository.
- `assets/annre-qr.png` / `.svg`: QR code for **https://github.com/Andy-B-123/AnnRE**.
- `assets/annre-qr.txt`: exact QR destination.

The repository must be published and publicly accessible before the QR is useful to your audience. Scan it on your phone after publishing and check it while signed out. Preserve the white border; use black on white, and allow ample slide space.

To regenerate the QR after installation:

```bash
python -m pip install '.[presentation]'
python scripts/make_qr.py https://github.com/Andy-B-123/AnnRE
```

For the local splash, open `docs/index.html` in a browser. A browser's full-page screenshot captures the full workflow. For an actual GitHub screenshot, publish first, open the repository at a readable zoom, hide personal browser/account details, and use your operating system's screenshot tool. The local splash is useful before publication but should not be described as a live repository screenshot.

Suggested spoken description:

> AnnRE—Annotation Refinement with Evidence—is an early prototype for using long-read RNA evidence to prioritise gene models for review. It complements transcript assembly with explicit checks of reference structures. Visual curation and tracked annotation changes are planned.

Suggested small-print attribution:

> Generated with OpenAI Codex and OpenAI models; reviewed and guided by Andreas Bachler.

Do not describe the 105 BSF candidates as validated corrections or claim an accuracy improvement over Bambu or ANNEXA. See [the BSF result summary](BSF_RESULTS.md) for a defensible comparison.
