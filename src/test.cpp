#include <bits/stdc++.h>
#include <sys/time.h>
using namespace std;

inline unsigned get_random() {
  static unsigned y = 2463534242;
  return y ^= (y ^= (y ^= y << 13) >> 17) << 5;
}

constexpr int ROW = 1 << 6;
constexpr int MAX_V = 1 << 9;
constexpr int MAX_KV = ROW * ROW;

int V, E, KV, KE, KR;
bool W[MAX_V][MAX_V];
int X[MAX_KV];
constexpr static int direction[] = {
    -ROW - 1, -ROW, -ROW + 1, -1, +1, +ROW - 1, +ROW, +ROW + 1,
};

void print() {
  for (int i = 1; i <= KR; ++i) {
    for (int j = 1; j <= KR; ++j) {
      int p = i * ROW + j;
      if (X[p] < V)
        fprintf(stderr, "%4d", X[i * ROW + j]);
      else
        fprintf(stderr, "    ");
    }
    fprintf(stderr, "\n");
  }
}

inline bool in(int p) {
  int r = p >> 6;
  int c = p & (ROW - 1);
  return 0 < r && r <= KR && 0 < c && c <= KR;
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
  {  // Hill Climbing
    int maxDistance = 5;
  start:
    for (int i = 0; i < MAX_KV; ++i) X[i] = V;
    vector<int> RemainVertex;
    for (int i = 1; i < V; ++i) RemainVertex.push_back(i);
    vector<int> used;
    {
      int t = KR / 2 * ROW + KR / 2;
      used.push_back(t);
      X[t] = 0;
    }
    vector<int> Position;
    while (RemainVertex.size()) {
      static int distance[MAX_KV][MAX_KV];
      memset(distance, 0x0f, sizeof(distance));
      for (int j = 0; j < MAX_KV; ++j) distance[j][j] = 0;
      Position.clear();
      for (int p : used) {
        Position.push_back(p);
        for (int d : direction)
          if (in(p + d) && X[p + d] == V) Position.push_back(p + d);
      }
      sort(Position.begin(), Position.end());
      Position.erase(unique(Position.begin(), Position.end()), Position.end());
      for (int p : Position) {
        if (X[p] < V) continue;
        for (int d : direction) {
          int n = p + d;
          if (in(n)) {
            distance[p][n] = 1;
            distance[n][p] = 1;
          }
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
      int center = -1;
      int value = 0;
      static int vertexDistance[MAX_V];
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
        for (int n : Position) {
          if (X[n] < V) continue;
          static int vd[MAX_V];
          memset(vd, 0x0f, sizeof(vd));
          for (int m : target) {
            int d = distance[n][m];
            if (vd[X[m]] > d) vd[X[m]] = d;
          }
          int t = 0;
          for (int m : targetVertex)
            if (vd[m] < maxDistance) t += 0xffff - vd[m];
          t = (t << 8) + (get_random() & ((1 << 8) - 1));
          if (value < t) {
            value = t;
            center = n;
            vertex = v;
            memcpy(vertexDistance, vd, sizeof(vd));
          }
        }
      }
      if (vertex == -1) {
        cerr << "bad" << endl;
        vertex = RemainVertex[0];
        for (int p : Position) {
          if (X[p] == V) {
            center = p;
            break;
          }
        }
        if (center == -1) {
          --maxDistance;
          goto start;
        }
      }
      // cerr << RemainVertex.size() << " " << used.size() << endl;
      // print();
      RemainVertex.erase(
          find(RemainVertex.begin(), RemainVertex.end(), vertex));
      X[center] = vertex;
      used.push_back(center);
      vector<int> tmp = used;
      for (int p : tmp) {
        if (W[vertex][X[p]] && vertexDistance[X[p]] < maxDistance &&
            vertexDistance[X[p]] == distance[center][p]) {
          int pos = p;
          while (pos != center) {
            for (int d : direction) {
              int n = pos + d;
              if (distance[center][pos] == distance[center][n] + 1 &&
                  (X[n] == vertex || X[n] == V)) {
                pos = n;
                if (X[pos] == V) {
                  X[pos] = vertex;
                  used.push_back(pos);
                }
              }
            }
          }
        }
        {
          bool empty = false;
          for (int d : direction) {
            empty |= X[p + d] == V;
          }
          if (!empty) used.erase(find(used.begin(), used.end(), p));
        }
      }
    }
  }
  {  // output
    vector<int> P[MAX_V];
    for (int i = 0; i < MAX_KV; ++i) {
      if (X[i] < V) P[X[i]].push_back((i / ROW - 1) * KR + i % ROW);
    }
    for (int i = 0; i < V; ++i) {
      printf("%d", (int)P[i].size());
      for (int v : P[i]) printf(" %d", v);
      printf("\n");
    }
  }
}
