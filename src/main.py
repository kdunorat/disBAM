from process_data import create_data
from del_search import DelSearch
import os
import time

if __name__ == '__main__':
    start = time.time()
    print('Loading Files...\n')
    header_log = 'Sample\tRegion\tSoft Clip\tSequence Size\tRegion Size\n'

    # Input files:
    absolute_path = os.path.dirname(__file__)
    files_list = list(filter(lambda s: '.bam' in s, os.listdir(f'{absolute_path}/input')))
    files_len = len(files_list)
    DelSearch.create_log(header_log, absolute_path)

    for (index, file) in enumerate(files_list):
        print(f'Searching deletions [{index+1}/{files_len}]')
        df, data, filename = create_data(file, absolute_path)

        run_del = DelSearch(df, data, filename)
        run_del.run_analysis(absolute_path)
    
    print(f"--- analysis ended in {(time.time() - start):.2f} seconds ---")
    print("You can access the log results on src/output/log.txt")
