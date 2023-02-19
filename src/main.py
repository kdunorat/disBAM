from process_data import create_data
from del_search import DelSearch
import time

if __name__ == '__main__':
    start = time.time()
    print('Loading Files...')

    # arquivo de entrada:
    file = "L144-247-232149084-NP.sorted.bam"
    df = create_data(file)

    run_del = DelSearch(df)
    print('\nSearching for deletions...')

    run_del.run_analysis()
    print(f"---{(time.time() - start)} seconds ---")
