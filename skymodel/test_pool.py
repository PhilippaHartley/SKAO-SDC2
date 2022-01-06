import multiprocessing as mp
import time
tstart = time.time()
def foo_pool(x):
    time.sleep(1)
    print (x*x)
    return x*x

result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)

def apply_async_with_callback():
    pool = mp.Pool(10)
    for i in range(100):
        print (i)
        pool.apply_async(foo_pool, args = (i, ), callback = log_result)
    pool.close()
    pool.join()
    print(result_list)

if __name__ == '__main__':
    apply_async_with_callback()
    tend = time.time()
    print (tstart-tend)