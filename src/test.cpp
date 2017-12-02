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

constexpr int ROW = 1 << 5;
constexpr int MAX_V = 1 << 9;
constexpr int MAX_KV = ROW * ROW;
constexpr int WALL = MAX_V - 1;

int V, E, KV, KE, KR;
uint8_t W[MAX_V][MAX_V];
uint8_t P[MAX_KV];
uint16_t X[MAX_KV];

int bestScore = 0;
uint16_t best[MAX_KV];

struct Node {
  int16_t pos;
  int16_t vertex;
  int16_t score;
  int value;
};
Node state[MAX_KV][MAX_V];

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
  int R = min(KR, ROW - 2);
  {  // Hill Climbing
    Node* bestVertex[MAX_KV];
    bool ok[MAX_KV];
    memset(ok, false, sizeof(ok));
    vector<int16_t> vertexes(MAX_V);
    vertexes.clear();
    vector<int16_t> positions(MAX_KV);
    positions.clear();
    for (int i = 0; i < MAX_KV; ++i) X[i] = WALL;
    for (int i = 0; i < V; ++i) vertexes.emplace_back(i);
    for (int i = 1; i <= R; ++i) {
      for (int j = 1; j <= R; ++j) {
        int p = i * ROW + j;
        ok[p] = true;
        X[p] = V;
        bestVertex[p] = &state[p][vertexes[0]];
        for (int k = 0; k < V; ++k) {
          state[p][k].pos = p;
          state[p][k].vertex = k;
          state[p][k].score = 0;
          state[p][k].value = INT_MIN;
        }
      }
    }
    int score = 0;
    auto update = [&](int n, int p) {
      if (X[p] != V) return;
      int value = -1;
      auto add = [&](Node& n, int v, int u) {
        if (u < V) {
          int t = W[v][u];
          if (t) {
            n.score += t;
            n.value += t << 14;
          } else {
            n.value -= 1 << 10;
          }
        }
      };
      for (int v : vertexes) {
        if (ok[p]) {
          state[p][v].score = 0;
          state[p][v].value = 0;
          add(state[p][v], v, X[p - ROW - 1]);
          add(state[p][v], v, X[p - ROW]);
          add(state[p][v], v, X[p - ROW + 1]);
          add(state[p][v], v, X[p - 1]);
          add(state[p][v], v, X[p + 1]);
          add(state[p][v], v, X[p + ROW - 1]);
          add(state[p][v], v, X[p + ROW]);
          add(state[p][v], v, X[p + ROW + 1]);
        } else {
          add(state[p][v], v, X[n]);
        }
        int t = state[p][v].value + (get_random() & ((1 << 10) - 1));
        if (value < t) {
          value = t;
          bestVertex[p] = &state[p][v];
        }
      }
      if (ok[p]) {
        ok[p] = false;
        positions.emplace_back(p);
      }
    };
    update(-1, (R / 2 + 1) * ROW + (R / 2 + 1));
    for (int i = 0; i < V; ++i) {
      Node* n = bestVertex[positions[0]];
      int value = INT_MIN;
      for (int p : positions) {
        if (bestVertex[p]->value < 0) {
          int value = INT_MIN;
          for (int v : vertexes) {
            if (value < state[p][v].value) {
              value = state[p][v].value;
              bestVertex[p] = &state[p][v];
            }
          }
        }
        int t = bestVertex[p]->value + (get_random() & ((1 << 10) - 1));
        if (value < t) {
          value = t;
          n = bestVertex[p];
        }
      }
      X[n->pos] = n->vertex;
      vertexes.erase(find(vertexes.begin(), vertexes.end(), n->vertex));
      positions.erase(find(positions.begin(), positions.end(), n->pos));
      score += n->score;
      for (int p : positions) state[p][n->vertex].value = -1000000;
      update(n->pos, n->pos - ROW - 1);
      update(n->pos, n->pos - ROW);
      update(n->pos, n->pos - ROW + 1);
      update(n->pos, n->pos - 1);
      update(n->pos, n->pos + 1);
      update(n->pos, n->pos + ROW - 1);
      update(n->pos, n->pos + ROW);
      update(n->pos, n->pos + ROW + 1);
    }
    if (bestScore < score) {
      bestScore = score;
      memcpy(best, X, sizeof(X));
    }
  }
  {  // output
    int P[MAX_V];
    for (int i = 0; i < MAX_KV; ++i) {
      if (best[i] < V) P[best[i]] = (i / ROW - 1) * KR + i % ROW;
    }
    for (int i = 0; i < V; ++i) {
      printf("1 %d\n", P[i]);
    }
  }
}