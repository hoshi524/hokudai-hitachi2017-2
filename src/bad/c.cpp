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
    if (V <= KR + 1) {
      for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
      --V;
      for (int i = 1; i <= V; ++i) {
        int p = ROW * i + 1, d = i & 1;
        for (int j = 1; j < V; ++j) {
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
      ++V;
      for (int i = 1; i < V; ++i) {
        X[i * ROW + V - 1] = V - 1;
      }
    } else {
      int size = 40, vertex;
      while (true) {
        for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
        vertex = 0;
        for (int i = 1; i + size <= KR + 1; i += size) {
          for (int j = 1; j <= KR; ++j) {
            if (size & 1) {
              if ((i + j) & 1) {
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
            } else {
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
        }
        if (vertex >= V) break;
        --size;
      }
      {
        for (int i = 0; i < MAX_KV; ++i) {
          int r = i >> 6;
          int c = i & (ROW - 1);
          if (r <= 1 || c <= 1 || r >= KR || c >= KR) continue;
          static int ROWA[] = {ROW, -ROW};
          for (int row : ROWA) {
            if (X[i] < vertex && X[i] == X[i + row + 1] &&
                X[i - row - 1] > vertex && X[i - row] < vertex &&
                X[i - row + 1] < vertex && X[i + row - 1] < vertex) {
              int D[] = {-row + 1, -row, -row - 1};
              for (int j = 0, p = i; j < size; ++j) {
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
              int D[] = {-row - 1, -row, -row + 1};
              for (int j = 0, p = i; j < size; ++j) {
                for (int d : D) {
                  if (X[p + d] > vertex) {
                    X[p + d] = X[p];
                    p += d;
                    break;
                  }
                }
              }
              i = 0;
            }
          }
        }
      }
      int16_t connect[MAX_V][MAX_V];
      auto connectSize = [&](int x) {
        int t;
        for (t = 0; connect[x][t] != -1; ++t) {
        }
        return t;
      };
      {  // connect
        auto calcConnectVertex = [&](int v) {
          static set<int> set;
          set.clear();
          for (int p = 0; p < MAX_KV; ++p) {
            if (X[p] == v) {
              for (int d : direction) {
                int n = p + d;
                if (X[n] < vertex && X[n] != v) {
                  set.insert(X[n]);
                }
              }
            }
          }
          int s = 0;
          for (int n : set) {
            connect[v][s++] = n;
          }
          return set.size();
        };
        auto calcConnect = [&]() {
          memset(connect, -1, sizeof(connect));
          int sum = 0;
          for (int v = 0; v < vertex; ++v) {
            sum += calcConnectVertex(v);
          }
          cerr << size << " " << (double)sum / vertex << endl;
        };
        calcConnect();
        {  // merge
          for (int t = V; t < vertex; ++t) {
            static int value[MAX_V][MAX_V];
            for (int i = 0; i < vertex; ++i) {
              for (int j = 0; j < vertex; ++j) {
                value[i][j] = INT_MIN;
              }
            }
            for (int r = 1; r <= KR; ++r) {
              for (int c = 1; c <= KR; ++c) {
                int p = r * ROW + c;
                if (X[p] > vertex) continue;
                for (int d : direction) {
                  int n = p + d;
                  if (X[n] > vertex) continue;
                  if (X[p] == X[n]) continue;
                  if (value[X[p]][X[n]] != INT_MIN) continue;
                  static set<int> set;
                  set.clear();
                  int ps, ns;
                  for (ps = 0; connect[X[p]][ps] != -1; ++ps) {
                    set.insert(connect[X[p]][ps]);
                  }
                  for (ns = 0; connect[X[n]][ns] != -1; ++ns) {
                    set.insert(connect[X[n]][ns]);
                  }
                  int v = -((ps + ns - set.size()) << 10) - (set.size() << 8) +
                          (get_random() & 0xff);
                  value[X[p]][X[n]] = v;
                  value[X[n]][X[p]] = v;
                }
              }
            }
            int v = INT_MIN, vmin, vmax;
            for (int i = 0; i < vertex; ++i) {
              for (int j = i + 1; j < vertex; ++j) {
                if (v < value[i][j]) {
                  v = value[i][j];
                  vmin = i;
                  vmax = j;
                }
              }
            }
            for (int i = 0; i < MAX_KV; ++i) {
              if (X[i] == vmax) X[i] = vmin;
            }
            calcConnectVertex(vmin);
          }
          int trans[MAX_KV];
          memset(trans, -1, sizeof(trans));
          for (int i = 0, v = 0; i < MAX_KV; ++i) {
            if (X[i] < vertex) {
              if (trans[X[i]] == -1) trans[X[i]] = v++;
              X[i] = trans[X[i]];
            }
          }
          vertex = V;
          calcConnect();
        }
        {
          if (X[ROW + KR] > vertex) {
            int x, v = 0xffff;
            for (int r = 1; r <= KR; ++r) {
              int p = r * ROW + KR - 1;
              if (X[p] > vertex) continue;
              int t = connectSize(X[p]);
              if (v > t) {
                v = t;
                x = X[p];
              }
            }
            for (int r = 1; r <= KR; ++r) {
              int p = r * ROW + KR - 1;
              if (X[p] > vertex) continue;
              X[p + 1] = x;
            }
            calcConnectVertex(x);
          }
          int r;
          for (r = 1; r <= KR; ++r)
            if (X[r * ROW + 1] > vertex) break;
          auto select = [&](int r, int &v1, int &v2) {
            int v = 0xffff;
            for (int c = 1; c <= KR; ++c) {
              int p = r * ROW + c;
              if (X[p] > vertex) continue;
              int t = connectSize(X[p]);
              if (v > t) {
                v = t;
                v1 = X[p];
              }
            }
            v = 0xffff;
            for (int c = 1; c <= KR; ++c) {
              int p = r * ROW + c;
              if (X[p] > vertex || v1 == X[p]) continue;
              int t = connectSize(X[p]);
              if (v > t) {
                v = t;
                v2 = X[p];
              }
            }
          };
          if (r + 1 <= KR) {
            int v1, v2;
            select(r - 1, v1, v2);
            for (int c = 1; c <= KR; ++c) {
              X[r * ROW + c] = c & 1 ? v1 : v2;
              X[(r + 1) * ROW + c] = c & 1 ? v2 : v1;
            }
          } else if (r <= KR) {
            int v1, v2;
            select(r - 1, v1, v2);
            for (int c = 1; c <= KR; ++c) {
              X[r * ROW + c] = v1;
            }
          }
          if (r + 3 <= KR) {
            for (int r = KR; r > 2; --r) {
              for (int c = 1; c <= KR; ++c) {
                X[r * ROW + c] = X[(r - 2) * ROW + c];
              }
            }
            int v1, v2;
            select(3, v1, v2);
            for (int c = 1; c <= KR; ++c) {
              X[ROW + c] = c & 1 ? v1 : v2;
              X[ROW + ROW + c] = c & 1 ? v2 : v1;
            }
          } else if (r + 2 <= KR) {
            for (int r = KR; r > 1; --r) {
              for (int c = 1; c <= KR; ++c) {
                X[r * ROW + c] = X[(r - 1) * ROW + c];
              }
            }
            int v1, v2;
            select(2, v1, v2);
            for (int c = 1; c <= KR; ++c) {
              X[ROW + c] = v1;
            }
          }
        }
        calcConnect();
      }
      int16_t x[MAX_V];
      int16_t rev[MAX_V];
      int16_t best[MAX_V];
      for (int i = 0; i < MAX_V; ++i) {
        x[i] = rev[i] = i < V ? i : MAX_V - 1;
      }
      auto value = [&](int v) {
        int16_t *c = connect[v];
        int t = 0;
        for (int i = 0; c[i] != -1; ++i) {
          t += W[x[v]][x[c[i]]];
        }
        return t;
      };
      int score = 0, bestScore = INT_MIN;
      if (V < 11) {
        do {
          score = 0;
          for (int i = 0; i < V; ++i) {
            score += value(i);
          }
          if (bestScore < score) {
            bestScore = score;
            memcpy(best, x, sizeof(x));
          }
        } while (next_permutation(x, x + V));
      } else {
        /**
         * 普通に性能悪い
         */
        constexpr int N = 5;
        constexpr double TIME_LIMIT = 10;
        while (TIME_LIMIT > timer.getElapsed()) {
          for (int k = 0; k < 0x1000; ++k) {
            static bool ok[MAX_V];
            memset(ok, true, sizeof(ok));
            static int B[N];
            for (int i = 0; i < N; ++i) {
              while (true) {
                int a = get_random() % vertex;
                if (ok[a]) {
                  ok[a] = false;
                  B[i] = a;
                  break;
                }
              }
            }
            static int T[N];
            static int C[N];
            for (int i = 0; i < N; ++i) {
              T[i] = x[B[i]];
              C[i] = i;
            }
            static int D[N];
            int score = INT_MIN;
            do {
              for (int i = 0; i < N; ++i) {
                x[B[i]] = T[C[i]];
              }
              int s = 0;
              for (int i = 0; i < N; ++i) {
                s += value(B[i]);
              }
              if (score < s) {
                score = s;
                memcpy(D, C, sizeof(C));
              }
            } while (next_permutation(C, C + N));
            for (int i = 0; i < N; ++i) {
              x[B[i]] = T[D[i]];
            }
          }
        }
        memcpy(best, x, sizeof(x));
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
