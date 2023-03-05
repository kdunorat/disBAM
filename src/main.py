from process_data import InputHandler, OutputHandler
from del_search import DelSearch
from depth_analyzer import DepthAnalyzer
import os
import time

if __name__ == '__main__':
    start = time.time()
    print('Loading Files...\n')

    # Input files:
    absolute_path = os.path.dirname(__file__)
    files_list = list(filter(lambda s: '.bam' in s, os.listdir(f'{absolute_path}/input')))
    files_len = len(files_list)
    header_log = 'Sample\tRegion\tSoft Clip\tSequence Size\tRegion Size\tCrossed\tReads Below 400\n'
    log = OutputHandler(absolute_path)
    log.create_log(header_log, 'w')

    for (index, file) in enumerate(files_list):
        print(f'Loading Data [{index + 1}/{files_len}]')
        input_data = InputHandler(file, absolute_path)
        df, data, filename, df_depth = input_data.create_data()

        # 1° feature
        print(f'Searching deletions [{index+1}/{files_len}]')
        run_del = DelSearch(df, data, filename)
        text = run_del.run_del_analysis()
        log.create_log(text, 'a')

        # 2° feature
        print(f'Analyzing depth [{index + 1}/{files_len}]')
        run_depth = DepthAnalyzer(df_depth, filename)
        text = run_depth.run_depth_analysis(7000, 400)
        log.create_log(text, 'a')

    print(f"--- analysis ended in {(time.time() - start):.2f} seconds ---")
    print("You can access the log results on src/output/log.txt")
