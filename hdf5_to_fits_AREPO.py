import argparse
import h5py
from pathlib import Path
import astropy.io.fits 

parser = argparse.ArgumentParser(description='Convert an HDF5 file to FITS')
parser.add_argument('--dataset', help='the dataset to read/write')
parser.add_argument('--output-path', type=Path, nargs='+',
	help='paths to write the outputs to')
parser.add_argument('path', type=Path, nargs='+',
	help='the paths to the snapshot files; any normal files are interpreted '
		'as single-file snapshots, while all files in any given directory '
		'are interpreted as a single multi-file snapshot')


def main():
	if not args.output_path:
		args.output_path = [None] * len(args.path)
	assert len(args.path) == len(args.output_path)
	# process paths
	for path_args in zip(args.path, args.output_path):
		process_path(*path_args)


def process_path(path, output_path):
	print(path)
	if not output_path:
		output_path = path.with_name(f'{path.stem}.fits')
	with h5py.File(path, 'r') as f:
		data = f[args.dataset][...]
	hdu = astropy.io.fits.PrimaryHDU(data.T)
	hdu.writeto(output_path)


if __name__ == '__main__':
	args = parser.parse_args()
	main()
