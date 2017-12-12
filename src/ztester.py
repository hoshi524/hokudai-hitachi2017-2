import subprocess
import queue
import threading


subprocess.call('mkdir out', shell=True)
subprocess.call('mkdir testcase', shell=True)
subprocess.call(
    'g++-5 -std=gnu++1y -O2 -o out/graph_generator.out toolkit/scripts/graph_generator.cpp', shell=True)
subprocess.call(
    'g++-5 -std=gnu++1y -O2 -o out/score_evaluator.out toolkit/scripts/score_evaluator.cpp', shell=True)
subprocess.call(
    'g++-5 -std=gnu++1y -O2 -o out/main.out src/main.cpp', shell=True)
subprocess.call(
    'g++-5 -std=gnu++1y -O2 -o out/test.out src/test.cpp', shell=True)
MAIN = './out/main.out'
TEST = './out/test.out'


def solve(command, seed):
    subprocess.call(
        '{0} < ./testcase/testcase{1} > ./testcase/result{1}'.format(command, seed), shell=True)
    score = int(subprocess.check_output(
        './out/score_evaluator.out ./testcase/testcase{0} ./testcase/result{0}'.format(seed), shell=True))
    return score


class State:
    main = 0
    test = 0
    lock = threading.Lock()

    def add(self, s, a, b):
        with self.lock:
            self.main += a
            self.test += b
            print('{}\t{}\t{}\t{}\t{}'.format(s, a, b, self.main, self.test))


scores = State()
q = queue.Queue()


def worker():
    while True:
        seed = q.get()
        if seed is None:
            break
        a = solve(MAIN, seed)
        b = solve(TEST, seed)
        scores.add(seed, a, b)
        q.task_done()


num_worker_threads = 3
threads = []
for i in range(num_worker_threads):
    t = threading.Thread(target=worker)
    t.start()
    threads.append(t)


N = 5000
for seed in range(0, N):
    subprocess.call(
        './out/graph_generator.out ./testcase/testcase{0} {0}'.format(seed), shell=True)
    q.put(seed)

# block until all tasks are done
q.join()

# stop workers
for i in range(num_worker_threads):
    q.put(None)
for t in threads:
    t.join()
