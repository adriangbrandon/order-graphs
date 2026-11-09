// Connection-aware greedy set cover on the De Bruijn graph (Section 6 of the paper), C++ version.
//
// Selects a set X of strings in [d]^d such that every requirement is covered and every string of X
// has a predecessor in X, then removes redundant strings (reverse delete) and writes X to a file.
// The order graph is obtained from X with decode_verify.py (proof of Lemma 13), which also checks
// completeness.
//
// Build:  g++ -O3 -march=native -std=c++17 -o greedy_ordergraph greedy_ordergraph.cpp
// Usage:  ./greedy_ordergraph d [--uni] [--seed S] [--alpha A] [--out file] [-v]
//
// Strings are encoded as base-d integers (first symbol most significant), symbols 0..d-1.
// A requirement (S, x) is encoded as x * 2^d + mask(S).

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>

using u64 = uint64_t;
using i64 = int64_t;

static int d;
static bool uni = false;
static u64 N;                  // d^d
static std::vector<u64> pw;    // powers of d

// Requirements covered by string w (at most 2d), appended to out.
static inline void reqs_of(u64 w, std::vector<u64> &out) {
    int digits[16];
    u64 t = w;
    for (int i = d - 1; i >= 0; --i) { digits[i] = int(t % d); t /= d; }
    int last = digits[d - 1];
    u64 seen = 0;
    for (int m = 1; m <= d; ++m) {
        int c = digits[d - m];
        if (seen >> c & 1) break;
        seen |= u64(1) << c;
        out.push_back(u64(last) << d | (seen & ~(u64(1) << last)));
        if (!uni && m > 1) out.push_back(u64(c) << d | (seen & ~(u64(1) << c)));
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s d [--uni] [--seed S] [--alpha A] [--out file] [-v]\n", argv[0]);
        return 1;
    }
    d = atoi(argv[1]);
    u64 seed = 0;
    double alpha = 1.0;
    std::string outname;
    bool verbose = false;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--uni") uni = true;
        else if (a == "--seed") seed = strtoull(argv[++i], nullptr, 10);
        else if (a == "--alpha") alpha = atof(argv[++i]);
        else if (a == "--out") outname = argv[++i];
        else if (a == "-v") verbose = true;
        else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 1; }
    }
    if (d < 2 || d > 9) { fprintf(stderr, "d must be in [2,9]\n"); return 1; }
    pw.assign(d + 1, 1);
    for (int i = 1; i <= d; ++i) pw[i] = pw[i - 1] * d;
    N = pw[d];
    if (outname.empty())
        outname = "strings_" + std::string(uni ? "uni" : "bi") + "_d" + std::to_string(d) + "_s" + std::to_string(seed) + ".txt";
    std::mt19937_64 rng(seed);
    auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

    // total requirements: d * 2^(d-1) (each requirement (S,x) with x not in S)
    u64 R = u64(d) << d;
    std::vector<char> covered(R, 0);
    u64 total = u64(d) << (d - 1), ncov = 0;

    std::vector<char> sel(N, 0);
    std::vector<u64> X;
    std::vector<std::vector<i64>> witness(d);
    for (int o = 0; o < d; ++o) witness[o].assign(pw[o], -1);

    std::vector<u64> buf;
    auto add = [&](u64 w) {
        if (sel[w]) return;
        sel[w] = 1;
        X.push_back(w);
        for (int o = 0; o < d; ++o) {
            i64 &slot = witness[o][w % pw[o]];
            if (slot < 0) slot = i64(w);
        }
        buf.clear();
        reqs_of(w, buf);
        for (u64 r : buf) if (!covered[r]) { covered[r] = 1; ++ncov; }
    };

    // initial ring: rotations of 0 1 ... d-1
    for (int r = 0; r < d; ++r) {
        u64 w = 0;
        for (int i = 0; i < d; ++i) w = w * d + u64((r + i) % d);
        add(w);
    }

    // candidates: strings with some uncovered requirement
    std::vector<u64> cand;
    cand.reserve(N);
    for (u64 s = 0; s < N; ++s) cand.push_back(s);

    std::vector<u64> stamp(R, 0);
    u64 stampv = 0;
    u64 path[16];
    int iter = 0;
    while (ncov < total) {
        ++iter;
        // filter candidates
        size_t k = 0;
        for (size_t i = 0; i < cand.size(); ++i) {
            u64 s = cand[i];
            if (sel[s]) continue;
            buf.clear();
            reqs_of(s, buf);
            bool useful = false;
            for (u64 r : buf) if (!covered[r]) { useful = true; break; }
            if (useful) cand[k++] = s;
        }
        cand.resize(k);

        double best = -1;
        u64 best_s = 0, ties = 0;
        for (u64 s : cand) {
            int o = d - 1;
            while (witness[o][s / pw[d - o]] < 0) --o;
            u64 w = u64(witness[o][s / pw[d - o]]);
            int len = 0;
            for (int j = o; j < d; ++j) {
                w = (w * d + (s / pw[d - 1 - j]) % d) % N;
                if (!sel[w]) path[len++] = w;
            }
            ++stampv;
            u64 gain = 0;
            for (int q = 0; q < len; ++q) {
                buf.clear();
                reqs_of(path[q], buf);
                for (u64 r : buf)
                    if (!covered[r] && stamp[r] != stampv) { stamp[r] = stampv; ++gain; }
            }
            double score = double(gain) / std::pow(double(len), alpha);
            if (score > best + 1e-12) { best = score; best_s = s; ties = 1; }
            else if (seed > 0 && score > best - 1e-12) {
                ++ties;
                if (rng() % ties == 0) best_s = s;
            }
        }
        // add the path of best_s
        {
            u64 s = best_s;
            int o = d - 1;
            while (witness[o][s / pw[d - o]] < 0) --o;
            u64 w = u64(witness[o][s / pw[d - o]]);
            for (int j = o; j < d; ++j) {
                w = (w * d + (s / pw[d - 1 - j]) % d) % N;
                add(w);
            }
        }
        if (verbose)
            fprintf(stderr, "  iter %d: %zu strings, %llu/%llu covered, %zu candidates, %.0fs\n", iter, X.size(),
                    (unsigned long long)ncov, (unsigned long long)total, cand.size(), elapsed());
    }
    size_t greedy_size = X.size();

    // reverse delete
    std::vector<u64> cnt(R, 0);
    for (u64 w : X) { buf.clear(); reqs_of(w, buf); for (u64 r : buf) ++cnt[r]; }
    auto has_other_pred = [&](u64 v, u64 except) {
        u64 suffix = v / d;  // v[1..d-1] as the last d-1 symbols of the predecessor
        for (int c = 0; c < d; ++c) {
            u64 p = u64(c) * pw[d - 1] + suffix;
            if (p != except && sel[p]) return true;
        }
        return false;
    };
    std::vector<u64> order = X;
    std::shuffle(order.begin(), order.end(), rng);
    bool changed = true;
    while (changed) {
        changed = false;
        for (u64 w : order) {
            if (!sel[w]) continue;
            buf.clear();
            reqs_of(w, buf);
            bool ok = true;
            for (u64 r : buf) if (cnt[r] < 2) { ok = false; break; }
            if (!ok) continue;
            for (int c = 0; c < d && ok; ++c) {
                u64 v = (w * d + c) % N;
                if (v != w && sel[v] && !has_other_pred(v, w)) ok = false;
            }
            if (!ok) continue;
            sel[w] = 0;
            for (u64 r : buf) --cnt[r];
            changed = true;
        }
    }
    std::vector<u64> final_set;
    for (u64 w : X) if (sel[w]) final_set.push_back(w);

    FILE *f = fopen(outname.c_str(), "w");
    fprintf(f, "# d=%d model=%s seed=%llu alpha=%g greedy=%zu pruned=%zu\n", d, uni ? "uni" : "bi",
            (unsigned long long)seed, alpha, greedy_size, final_set.size());
    for (u64 w : final_set) {
        char s[16];
        u64 t = w;
        for (int i = d - 1; i >= 0; --i) { s[i] = char('1' + t % d); t /= d; }
        s[d] = 0;
        fprintf(f, "%s\n", s);
    }
    fclose(f);
    printf("d=%d %s seed=%llu: greedy %zu -> pruned %zu strings, %.1fs -> %s\n", d, uni ? "uni" : "bi",
           (unsigned long long)seed, greedy_size, final_set.size(), elapsed(), outname.c_str());
    return 0;
}
