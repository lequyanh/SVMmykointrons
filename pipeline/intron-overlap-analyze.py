import csv
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    filename='intron-overlap-analyze.log',
    filemode='w'
)


def process_with_label():
    for scaffold, positions in scaffolds_dict.items():
        last_end, last_start, last_label = 0, 0, 0

        for start, end, label, _ in positions:
            start = int(start)
            end = int(end)

            if start < last_end:
                overlap_ratio = (last_end - start) / (last_end - last_start)
                logging.info(
                    f'{scaffold}-{last_start}-{last_end}---{start}-{end}-{overlap_ratio}--={last_label};{label}')

            last_end = end
            last_start = start
            last_label = label


def process_without_label():
    for scaffold, positions in scaffolds_dict.items():
        last_end, last_start = 0, 0

        for start, end in positions:
            start = int(start)
            end = int(end)

            if start < last_end:
                overlap_ratio = (last_end - start) / (last_end - last_start)
                logging.info(f'{scaffold}-{last_start}-{last_end}---{start}-{end}-{overlap_ratio}')

            last_end = end
            last_start = start


if __name__ == "__main__":
    filename = sys.argv[1]  # intron-result.csv
    logging.info(f'Analyzing intron suspects (the pipeline output) in file {filename}')

    with open(filename, 'r') as f:
        content = csv.reader(f, delimiter=';')

        scaffolds_dict = dict()
        for i, row in enumerate(content):
            if i == 0:
                has_label = len(row) > 3  # scaffold;position;label;prediction
                continue

            scaff = row[0]
            positions = scaffolds_dict.get(scaff, [])
            positions.append((row[1:]))

            scaffolds_dict[scaff] = positions

    if has_label:
        process_with_label()
    else:
        process_without_label()
