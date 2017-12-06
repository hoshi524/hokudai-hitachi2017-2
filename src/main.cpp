#include <bits/stdc++.h>
#include <sys/time.h>
using namespace std;

class Timer {
 public:
  void restart();
  double getElapsed();

  Timer();

 private:
  static double rdtsc_per_sec_inv;

  double getTimeOfDay();
  unsigned long long int getCycle();

  double start_time;
  unsigned long long int start_clock;
};
double Timer::rdtsc_per_sec_inv = -1;

inline double Timer::getElapsed() {
  if (rdtsc_per_sec_inv != -1)
    return (double)(getCycle() - start_clock) * rdtsc_per_sec_inv;

  const double RDTSC_MEASUREMENT_INTERVAL = 0.1;
  double res = getTimeOfDay() - start_time;
  if (res <= RDTSC_MEASUREMENT_INTERVAL) return res;

  rdtsc_per_sec_inv = 1.0 / (getCycle() - start_clock);
  rdtsc_per_sec_inv *= getTimeOfDay() - start_time;
  return getElapsed();
}

inline void Timer::restart() {
  start_time = getTimeOfDay();
  start_clock = getCycle();
}

Timer::Timer() { restart(); }

inline double Timer::getTimeOfDay() {
  timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + tv.tv_usec * 0.000001;
}

inline unsigned long long int Timer::getCycle() {
  unsigned int low, high;
  __asm__ volatile("rdtsc" : "=a"(low), "=d"(high));
  return ((unsigned long long int)low) | ((unsigned long long int)high << 32);
}

Timer timer;

inline unsigned get_random() {
  static unsigned y = 2463534242;
  return y ^= (y ^= (y ^= y << 13) >> 17) << 5;
}

constexpr int ROW = 1 << 6;
constexpr int MAX_V = 1 << 10;
constexpr int MAX_KV = ROW * ROW;

int V, E, KV, KE, KR;
int8_t W[MAX_V][MAX_V];
int X[MAX_KV];
constexpr static int direction[] = {
    -ROW - 1, -ROW, -ROW + 1, -1, +1, +ROW - 1, +ROW, +ROW + 1,
};

void print() {
  for (int i = 1; i <= KR; ++i) {
    for (int j = 1; j <= KR; ++j) {
      int p = i * ROW + j;
      if (X[p] < MAX_V - 1)
        fprintf(stderr, "%4d", X[i * ROW + j]);
      else
        fprintf(stderr, "    ");
    }
    fprintf(stderr, "\n");
  }
}

int main() {
  {  // input
    memset(W, 0, sizeof(W));
    scanf("%d%d\n", &V, &E);
    int u, v;
    for (int i = 0; i < E; ++i) {
      scanf("%d%d\n", &u, &v);
      --u;
      --v;
      W[u][v] = 1;
      W[v][u] = 1;
    }
    scanf("%d%d\n", &KV, &KE);
    KR = sqrt(KV);
  }
  {  // Annealing
    if (V <= KR) {
      for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
      for (int i = 1; i <= V; ++i) {
        int p = ROW * i + 1, d = i & 1;
        for (int j = 1; j <= V; ++j) {
          X[p] = i - 1;
          int n;
          if (d) {
            n = p + 1 + ROW;
          } else {
            n = p + 1 - ROW;
          }
          if (n < ROW || n >= (V + 1) * ROW) {
            n = p + 1;
            d ^= 1;
          } else if (X[n] != MAX_V - 1) {
            n = p + 1;
          }
          p = n;
        }
      }
    } else {
      int size = 40, vertex;
      while (true) {
        for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
        vertex = 0;
        for (int i = 1; i + size <= KR + 1; i += size) {
          for (int j = 1; j <= KR; ++j) {
            if (j & 1) {
              if (j + size <= KR + 1) {
                for (int k = 0, p = i * ROW + j; k < size; ++k) {
                  X[p] = vertex;
                  p += ROW + 1;
                }
                ++vertex;
              }
            } else {
              if (j >= size) {
                for (int k = 0, p = i * ROW + j; k < size; ++k) {
                  X[p] = vertex;
                  p += ROW - 1;
                }
                ++vertex;
              }
            }
          }
        }
        if (vertex >= V) break;
        --size;
      }
      {
        bool ok = true;
        while (ok) {
          ok = false;
          for (int i = 0; i < MAX_KV; ++i) {
            int r = i >> 6;
            int c = i & (ROW - 1);
            if (r < 2 || c < 2 || r >= KR - 1 || c >= KR - 1) continue;
            static int ROWA[] = {ROW, -ROW};
            for (int row : ROWA) {
              if (X[i] < vertex && X[i] == X[i + row + 1] &&
                  X[i - row - 1] > vertex && X[i - row] < vertex &&
                  X[i - row + 1] < vertex && X[i + row - 1] < vertex) {
                ok = true;
                for (int j = 0, p = i; j < size; ++j) {
                  int D[] = {-row + 1, -row, -row - 1};
                  for (int d : D) {
                    if (X[p + d] > vertex) {
                      X[p + d] = X[p];
                      p += d;
                      break;
                    }
                  }
                }
              }
              if (X[i] < vertex && X[i] == X[i + row - 1] &&
                  X[i - row + 1] > vertex && X[i - row] < vertex &&
                  X[i - row - 1] < vertex && X[i + row + 1] < vertex) {
                ok = true;
                for (int j = 0, p = i; j < size; ++j) {
                  int D[] = {-row - 1, -row, -row + 1};
                  for (int d : D) {
                    if (X[p + d] > vertex) {
                      X[p + d] = X[p];
                      p += d;
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
      int16_t connect[MAX_V][MAX_V];
      memset(connect, -1, sizeof(connect));
      int sum = 0;
      for (int v = 0; v < vertex; ++v) {
        static set<int> set;
        set.clear();
        for (int p = 0; p < MAX_KV; ++p) {
          if (X[p] == v) {
            for (int d : direction) {
              int n = p + d;
              if (X[n] != V && X[n] != v) {
                set.insert(X[n]);
              }
            }
          }
        }
        int s = 0;
        for (int n : set) {
          connect[v][s++] = n;
        }
        sum += set.size();
      }
      cerr << size << " " << (double)sum / vertex << endl;
      int16_t x[MAX_V];
      int16_t best[MAX_V];
      for (int i = 0; i < MAX_V; ++i) {
        x[i] = i < V ? i : MAX_V - 1;
      }
      constexpr double TIME_LIMIT = 1.9;
      constexpr int LOG_SIZE = 1 << 10;
      double log_d[LOG_SIZE];
      uint8_t log_[LOG_SIZE];
      for (int i = 0; i < LOG_SIZE; ++i) {
        log_d[i] = -0.5 * log((i + 0.5) / LOG_SIZE) / TIME_LIMIT;
      }
      int score = 0, bestScore = 0;
      while (true) {
        double time = TIME_LIMIT - timer.getElapsed();
        if (time < 0) break;
        for (int i = 0; i < LOG_SIZE; ++i)
          log_[i] = min(10.0, round(log_d[i] * time));
        for (int t = 0; t < 0x10000; ++t) {
          int a = get_random() % vertex;
          int b = get_random() % vertex;
          if (a == b) continue;
          auto value = [&](int v) {
            int t = 0;
            for (int i = 0; connect[v][i] != -1; ++i) {
              t += W[x[v]][x[connect[v][i]]];
            }
            return t;
          };
          int pv = value(a) + value(b);
          swap(x[a], x[b]);
          int nv = value(a) + value(b);
          int d = pv - nv;
          if (d > log_[get_random() & (LOG_SIZE - 1)]) {
            swap(x[a], x[b]);
          } else {
            score -= d;
            if (bestScore < score) {
              bestScore = score;
              memcpy(best, x, sizeof(x));
            }
          }
        }
      }
      int T[MAX_KV];
      memcpy(T, X, sizeof(T));
      for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
      for (int i = 0; i < MAX_V; ++i) {
        for (int j = 0; j < MAX_KV; ++j) {
          if (T[j] == i) {
            X[j] = best[i];
          }
        }
      }
    }
  }
  {  // output
    vector<int> P[MAX_V];
    for (int i = 0; i < MAX_KV; ++i) {
      if (X[i] < MAX_V - 1) P[X[i]].push_back((i / ROW - 1) * KR + i % ROW);
    }
    for (int i = 0; i < V; ++i) {
      printf("%d", (int)P[i].size());
      for (int v : P[i]) printf(" %d", v);
      printf("\n");
    }
  }
}
