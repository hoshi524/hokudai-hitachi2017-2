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
constexpr int MAX_V = 1 << 9;
constexpr int MAX_KV = ROW * ROW;

int V, E, KV, KE, KR;
int W[MAX_V][MAX_V];
int X[MAX_KV];
constexpr static int direction[] = {
    -ROW - 1, -ROW, -ROW + 1, -1, +1, +ROW - 1, +ROW, +ROW + 1,
};

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
  {  // Hill Climbing
  start:
    for (int i = 0; i < MAX_KV; ++i) X[i] = V;
    vector<int> RemainVertex;
    for (int i = 0; i < V; ++i) RemainVertex.push_back(i);
    vector<int> Position;
    Position.push_back(KR / 2 * ROW + KR / 2);
    vector<int> used;
    for (int i = 0; i < V; ++i) {
      static int distance[MAX_KV][MAX_KV];
      memset(distance, 0x0f, sizeof(distance));
      for (int j = 0; j < V; ++j) distance[j][j] = 0;
      for (int p : Position) {
        for (int d : direction) {
          int n = p + d;
          if (X[n] < V) continue;
          distance[p][n] = 1;
          distance[n][p] = 1;
        }
      }
      for (int k : Position) {
        if (X[k] < V) continue;
        for (int i : Position) {
          for (int j : Position) {
            int d = distance[i][k] + distance[k][j];
            if (distance[i][j] > d) distance[i][j] = d;
          }
        }
      }
      int vertex = -1;
      int value = 0;
      for (int v : RemainVertex) {
        static vector<int> target;
        static set<int> targetVertex;
        target.clear();
        targetVertex.clear();
        for (int p : used) {
          if (W[v][X[p]]) {
            target.push_back(p);
            targetVertex.insert(X[p]);
          }
        }
        if (target.empty()) continue;
        for (int n : target) {
          static int VertexDistance[MAX_V];
          memset(VertexDistance, 0x0f, sizeof(VertexDistance));
          for (int m : target) {
            int d = distance[n][m];
            if (VertexDistance[X[m]] > d) VertexDistance[X[m]] = d;
          }
          int v = 0;
          for (int n : targetVertex)
            if (VertexDistance[n] < 0xff) v += 0xffff - VertexDistance[n];
          if (value < v) {
            value = v;
            vertex = v;
          }
        }
      }
      if (vertex == -1) goto start;
    }
  }
}
