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
constexpr int MAX_V = 1 << 11;
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
    if (V == 2) {
      for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
      X[ROW + 1] = 0;
      X[ROW + 2] = 1;
    } else if (V <= KR + 1) {
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
      int size, vertex;
      static int16_t connect[MAX_V][MAX_V];
      static int16_t CS[MAX_V];
      auto calcConnectVertex = [&](int v) {
        static vector<int> set(MAX_V);
        set.clear();
        for (int p = 0; p < MAX_KV; ++p) {
          if (X[p] == v) {
            for (int d : direction) {
              int n = p + d;
              if (X[n] < vertex && X[n] != v) set.emplace_back(X[n]);
            }
          }
        }
        sort(set.begin(), set.end());
        set.erase(unique(set.begin(), set.end()), set.end());
        memset(connect[v], -1, sizeof(connect[v]));
        int s = 0;
        for (int n : set) {
          connect[v][s++] = n;
        }
        CS[v] = s;
        return s;
      };
      auto calcConnect = [&]() {
        int sum = 0;
        for (int v = 0; v < vertex; ++v) {
          sum += calcConnectVertex(v);
        }
        // cerr << size << " " << (double)sum / vertex << endl;
        return sum;
      };
      int CX[MAX_KV], cscore = 0;
      for (int q = 0; q < 5; ++q) {
        size = 40;
        vertex = 0;
        auto setX = [&](int p, int t) {
          for (int k = 0; k < size; ++k) {
            X[p] = vertex;
            p += ROW + t;
          }
          ++vertex;
        };
        while (true) {
          for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
          vertex = q;
          for (int i = 1 + min(2, q); i + size <= KR + 1 - max(0, q - 2);
               i += size) {
            for (int j = 1; j <= KR; ++j) {
              int p = i * ROW + j;
              if (size & 1) {
                if ((i + j) & 1) {
                  if (j + size <= KR + 1) setX(p, 1);
                } else {
                  if (j >= size) setX(p, -1);
                }
              } else {
                if (j & 1) {
                  if (j + size <= KR + 1) setX(p, 1);
                } else {
                  if (j >= size) setX(p, -1);
                }
              }
            }
          }
          if (vertex >= V) break;
          --size;
          if (size == 0) goto XSetEnd;
        }
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
        if (X[(1 + +min(2, q)) * ROW + KR] > vertex) {
          for (int r = 1; r <= KR; ++r) {
            if (X[r * ROW + KR - 1] < vertex) X[r * ROW + KR] = vertex;
          }
          ++vertex;
        }
        if (q == 1) {
          for (int i = 1; i <= KR; ++i) {
            X[1 * ROW + i] = 0;
          }
        } else if (q == 2) {
          for (int i = 1; i <= KR; ++i) {
            X[1 * ROW + i] = (i & 1) ? 0 : 1;
            X[2 * ROW + i] = (i & 1) ? 1 : 0;
          }
        } else if (q == 3) {
          for (int i = 1; i <= KR; ++i) {
            X[1 * ROW + i] = (i & 1) ? 0 : 1;
            X[2 * ROW + i] = (i & 1) ? 1 : 0;
            X[(KR - 0) * ROW + i] = 2;
          }
        } else if (q == 4) {
          for (int i = 1; i <= KR; ++i) {
            X[1 * ROW + i] = (i & 1) ? 0 : 1;
            X[2 * ROW + i] = (i & 1) ? 1 : 0;
            X[(KR - 1) * ROW + i] = (i & 1) ? 2 : 3;
            X[(KR - 0) * ROW + i] = (i & 1) ? 3 : 2;
          }
        }
        for (int r = 2; r <= KR; ++r) {
          for (int c = 1; c <= KR; ++c) {
            int p = r * ROW + c;
            if (X[p] < vertex) continue;
            if (X[p - ROW + 1] < vertex &&
                (X[p - ROW + 1] == X[p - ROW - ROW + 2] ||
                 X[p - ROW + 1] == X[p - ROW - ROW + 1] ||
                 X[p - ROW + 1] == X[p - ROW - ROW + 0])) {
              X[p] = X[p - ROW + 1];
            }
            if (X[p - ROW - 1] < vertex &&
                (X[p - ROW - 1] == X[p - ROW - ROW - 2] ||
                 X[p - ROW - 1] == X[p - ROW - ROW - 1] ||
                 X[p - ROW + 1] == X[p - ROW - ROW - 0])) {
              X[p] = X[p - ROW - 1];
            }
            if (X[p] > vertex) X[p] = X[p - ROW];
          }
        }
        calcConnect();
        {  // merge
          for (int t = V; t < vertex; ++t) {
            static bool used[MAX_V][MAX_V];
            memset(used, false, sizeof(used));
            int v = INT_MIN, vmin, vmax;
            for (int r = 1; r <= KR; ++r) {
              for (int c = 1; c <= KR; ++c) {
                int p = r * ROW + c;
                if (X[p] > vertex) continue;
                for (int d : direction) {
                  int n = p + d;
                  if (X[n] > vertex) continue;
                  if (X[p] == X[n]) continue;
                  if (used[X[p]][X[n]]) continue;
                  used[X[p]][X[n]] = true;
                  used[X[n]][X[p]] = true;
                  static vector<int> set(MAX_V);
                  set.clear();
                  int ps, ns;
                  for (ps = 0; connect[X[p]][ps] != -1; ++ps) {
                    set.emplace_back(connect[X[p]][ps]);
                  }
                  for (ns = 0; connect[X[n]][ns] != -1; ++ns) {
                    set.emplace_back(connect[X[n]][ns]);
                  }
                  sort(set.begin(), set.end());
                  int s = unique(set.begin(), set.end()) - set.begin();
                  int g =
                      -((ps + ns - s) << 16) - (s << 8) + (get_random() & 0xff);
                  if (v < g) {
                    v = g;
                    vmin = min(X[p], X[n]);
                    vmax = max(X[p], X[n]);
                  }
                }
              }
            }
            for (int i = 0; i < MAX_KV; ++i) {
              if (X[i] == vmax) X[i] = vmin;
            }
            for (int i = 0; i < vertex; ++i) {
              bool a = false, b = false;
              for (int j = 0; connect[i][j] != -1; ++j) {
                a |= connect[i][j] == vmax;
                b |= connect[i][j] == vmin;
              }
              if (a) {
                for (int j = 0; connect[i][j] != -1; ++j) {
                  if (connect[i][j] == vmax) {
                    if (b) {
                      connect[i][j] = connect[i][--CS[i]];
                      connect[i][CS[i]] = -1;
                    } else {
                      connect[i][j] = vmin;
                    }
                  }
                }
              }
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
        }
        int s = calcConnect();
        if (cscore < s) {
          cscore = s;
          memcpy(CX, X, sizeof(X));
        }
        if (size == 1) goto XSetEnd;
      }
    XSetEnd:
      memcpy(X, CX, sizeof(X));
      calcConnect();
      static int16_t x[MAX_V];
      static int16_t rev[MAX_V];
      static int16_t best[MAX_V];
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
      if (V < 12) {
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
        static int16_t WA[MAX_V][MAX_V];
        static int16_t WS[MAX_V];
        static int16_t P[MAX_V];
        memset(WS, 0, sizeof(WS));
        for (int i = 0; i < vertex; ++i) {
          P[i] = value(i);
          for (int j = 0; j < vertex; ++j) {
            if (W[i][j]) WA[i][WS[i]++] = j;
          }
        }
        constexpr double TIME_LIMIT = 2.9;
        constexpr int LOG_SIZE = 1 << 10;
        static double log_d[LOG_SIZE];
        static uint8_t log_[LOG_SIZE];
        for (int i = 0; i < LOG_SIZE; ++i) {
          log_d[i] = -0.5 * log((i + 0.5) / LOG_SIZE) / TIME_LIMIT;
        }
        while (true) {
          double time = TIME_LIMIT - timer.getElapsed();
          if (time < 0) break;
          for (int i = 0; i < LOG_SIZE; ++i)
            log_[i] = min(10.0, round(log_d[i] * time));
          for (int t = 0; t < 0x100; ++t) {
            for (int a = 0; a < vertex; ++a) {
              int z = rev[WA[x[a]][get_random() % WS[x[a]]]];
              int b = connect[z][get_random() % CS[z]];
              if (a == b) continue;
              int pv = P[a] + P[b];
              swap(x[a], x[b]);
              int va = value(a);
              int vb = value(b);
              int d = pv - va - vb;
              if (d > log_[get_random() & (LOG_SIZE - 1)]) {
                swap(x[a], x[b]);
              } else {
                rev[x[a]] = a;
                rev[x[b]] = b;
                auto diff = [&](int p, int v) {
                  int16_t *c = connect[v];
                  for (int i = 0; c[i] != -1; ++i)
                    P[c[i]] += W[x[v]][x[c[i]]] - W[p][x[c[i]]];
                };
                diff(x[a], b);
                diff(x[b], a);
                P[a] = va;
                P[b] = vb;
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
      static int T[MAX_KV];
      memcpy(T, X, sizeof(T));
      for (int i = 0; i < MAX_KV; ++i) X[i] = MAX_V - 1;
      for (int i = 0; i < vertex; ++i) {
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