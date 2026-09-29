#!/usr/bin/env python3
"""Generate a conference QR with its exact destination recorded alongside it."""
import argparse
from pathlib import Path
from urllib.parse import urlparse
import qrcode
import qrcode.image.svg

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('url', help='Public HTTPS page to open, not a .git clone URL')
p.add_argument('--out', type=Path, default=Path('docs/assets/annre-qr'))
a = p.parse_args()
u = urlparse(a.url)
if u.scheme != 'https' or not u.netloc or u.path.endswith('.git'):
    p.error('Use a public HTTPS page URL without the .git suffix.')
a.out.parent.mkdir(parents=True, exist_ok=True)
qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_M, box_size=16, border=4)
qr.add_data(a.url)
qr.make(fit=True)
qr.make_image(fill_color='black', back_color='white').save(str(a.out)+'.png')
qr.make_image(image_factory=qrcode.image.svg.SvgPathImage).save(str(a.out)+'.svg')
Path(str(a.out)+'.txt').write_text(a.url+'\n')
print(f'QR destination: {a.url}')
