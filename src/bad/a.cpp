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
        if (size & 1) {
          for (int r = 1; r < KR; ++r) {
            for (int c = 1; c < KR; ++c) {
              int i = r * ROW + c;
              static int K[] = {1, -1};
              for (int k : K) {
                if (X[i] < vertex && X[i - k] > vertex && X[i + k] < vertex &&
                    X[i + ROW] < vertex && X[i - ROW] < vertex &&
                    X[i - ROW - k] < vertex && X[i + ROW - k] < vertex &&
                    X[i] != X[i - ROW - k] && X[i] != X[i + ROW - k]) {
                  X[i - k] = X[i];
                  i -= k;
                  static int ROWA[] = {ROW, -ROW};
                  for (int row : ROWA) {
                    int di = 3;
                    int D[] = {-row + k, -row, -row - k};
                    for (int j = 0, p = i; j < size; ++j) {
                      for (int j = 0; j < di; ++j) {
                        int d = D[j];
                        if (X[p + d] > vertex &&
                            (X[p + d + 1] < vertex || X[p + d - 1] < vertex) &&
                            (X[p + d + row] < vertex ||
                             X[p + d - row] < vertex)) {
                          X[p + d] = X[p];
                          p += d;
                          di = j + 1;
                          break;
                        }
                      }
                    }
                  }
                  r = 1;
                  c = 1;
                }
              }
            }
          }
          // 適当に埋める
          for (int r = 1; r <= KR; ++r) {
            for (int c = 1; c <= KR; ++c) {
              int p = r * ROW + c;
              if (X[p] > vertex) {
                int t = 0;
                if (X[p - 1] < vertex) ++t;
                if (X[p + 1] < vertex) ++t;
                if (X[p - ROW] < vertex) ++t;
                if (X[p + ROW] < vertex) ++t;
                if (t >= 2) {
                  int v = 0;
                  static int D[] = {ROW - 1, ROW + 1, -ROW - 1, -ROW + 1};
                  for (int d : D) {
                    int n = p + d;
                    if (X[n] < vertex) {
                      int a = get_random();
                      if (v < a) {
                        v = a;
                        t = n;
                      }
                    }
                  }
                  X[p] = X[t];
                  r = 1;
                  c = 1;
                }
              }
            }
          }
        } else {
          for (int i = 0; i < MAX_KV; ++i) {
            int r = i >> 6;
            int c = i & (ROW - 1);
            if (r < 2 || c < 2 || r >= KR - 1 || c >= KR - 1) continue;
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
          calcConnect();
        }
      }
      int16_t x[MAX_V];
      int16_t best[MAX_V];
      for (int i = 0; i < MAX_V; ++i) {
        x[i] = i < V ? i : MAX_V - 1;
      }
      auto value = [&](int v) {
        int t = 0;
        for (int i = 0; connect[v][i] != -1; ++i) {
          t += W[x[v]][x[connect[v][i]]];
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
        int sorted[MAX_V];
        int S[MAX_V];
        int T[MAX_V];
        for (int i = 0; i < vertex; ++i) {
          sorted[i] = i;
          S[i] = connectSize(i);
        }
        sort(sorted, sorted + vertex,
             [S](int &a, int &b) { return S[a] > S[b]; });
        int range[MAX_V][2];
        bool perfect = true;
        {
          auto calcT = [&]() {
            for (int i = 0; i < vertex; ++i) {
              int t = 0;
              for (int j = 0; j < vertex; ++j) t += W[i][j];
              T[i] = t;
            }
          };
          calcT();
          sort(T, T + vertex, [](int &a, int &b) { return a > b; });
          for (int i = 0; i < vertex; ++i) {
            perfect &= S[sorted[i]] >= T[i];
          }
          calcT();
          if (perfect) {
            /**
             * 終域の最小次数 > 定義域の最大次数
             * でも完全グラフにほぼならない
             * から、そこで遷移先を狭めようとしてもあまり効果がなさそうだった
             */
            bool left[MAX_V];
            memset(left, true, sizeof(left));
            for (int i = 0; i < vertex; ++i) {
              int j;
              for (j = 0; j < vertex; ++j) {
                if (left[j]) break;
              }
              range[i][0] = j;
              for (j = 0; j < vertex; ++j) {
                if (S[sorted[j]] < T[i]) break;
              }
              range[i][1] = j;
              for (j = j - 1; !left[j]; --j) {
              }
              left[j] = false;
            }
            {
              for (int i = 0; i < vertex; ++i) {
                cerr << range[i][0] << "~" << range[i][1] << "\t";
              }
              cerr << endl;
            }
          } else {
            for (int i = 0; i < vertex; ++i) {
              range[i][0] = 0;
              range[i][1] = vertex;
            }
          }
        }
        for (int i = 0; i < vertex; ++i) {
          range[i][1] -= range[i][0];
        }
        constexpr double TIME_LIMIT = 1.9;
        constexpr int LOG_SIZE = 1 << 10;
        double log_d[LOG_SIZE];
        uint8_t log_[LOG_SIZE];
        for (int i = 0; i < LOG_SIZE; ++i) {
          log_d[i] = -0.5 * log((i + 0.5) / LOG_SIZE) / TIME_LIMIT;
        }
        while (true) {
          double time = TIME_LIMIT - timer.getElapsed();
          if (time < 0) break;
          for (int i = 0; i < LOG_SIZE; ++i)
            log_[i] = min(10.0, round(log_d[i] * time));
          for (int t = 0; t < 0x10; ++t) {
            for (int a = 0; a < vertex; ++a) {
              int b = sorted[range[x[a]][0] + get_random() % range[x[a]][1]];
              if (a == b) continue;
              assert(!(perfect && S[b] < T[x[a]]));
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
