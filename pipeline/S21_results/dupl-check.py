import csv
import logging

logging.basicConfig(
    level=logging.INFO,
    filename='dupl-check.log',
    filemode='w'
)

if __name__ == "__main__":
    filename = 'intron-result.csv'

    with open(filename, 'r') as f:
        content = csv.reader(f, delimiter=';')

        scaffolds_dict = dict()
        for i, row in enumerate(content):
            if i == 0:
                continue

            scaff = row[0]
            positions = scaffolds_dict.get(scaff, [])
            positions.append((row[1], row[2]))

            scaffolds_dict[scaff] = positions

    for scaffold, positions in scaffolds_dict.items():
        last_end = 0
        last_start = 0

        for start, end in positions:
            start = int(start)
            end = int(end)

            if start < last_end:
                overlap_ratio = (last_end - start)/(last_end - last_start)
                logging.info(f'{scaffold}-{last_start}-{last_end}---{start}-{end}-{overlap_ratio}')

            last_end = end
            last_start = start
