#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <iomanip>

using std::vector;
using std::pair;

using T = int;
using MT = long long;

class T_task {
    vector<vector<T>> t_body;
    vector<vector<T>> t_help_body;
    vector<vector<T>> coasts;
    vector<pair<size_t, size_t>> find_base_sln(vector<T> producer, vector<T> consumer) {
        // Методом верхнего левого квадрата
        vector<pair<size_t, size_t>> basis_elems;
        T ostatok = producer[1];
        for (size_t i = 1, j = 1; true;) {
            if (ostatok - consumer[i] > 0) {
                t_body[i][j] = consumer[i];
                basis_elems.emplace_back(i, j);
                ostatok -= consumer[i];
                ++i;
            }
            else if (ostatok - consumer[i] < 0) {
                t_body[i][j] = ostatok;
                basis_elems.emplace_back(i, j);
                consumer[i] -= ostatok;
                ++j;
                ostatok = producer[j];
            }
            else {//ostatok == 0
                t_body[i][j] = ostatok;
                basis_elems.emplace_back(i, j);
                ++j;
                ++i;
                if (j >= producer.size()) break;
                ostatok = producer[j];
            }
        }
        return basis_elems;
    }
    
    void append_basis(vector<pair<size_t, size_t>>& basis) {
        size_t size = basis.size();
        for (size_t i = 0; i != size - 1; ++i) {
            if (not((basis[i].first == basis[i + 1].first) or (basis[i].second == basis[i + 1].second))) {
                basis.emplace_back(basis[i].first, basis[i].second + 1);
            }
        }
        sort(basis);
    }
        
    void find_potential(vector<MT>& alpha, vector<MT>& beta, vector<pair<size_t, size_t>> basis) {
        alpha[1] = coasts[1][1];
        beta[1] = 0;
        
        size_t a = 1, b = 1;

        for (size_t i = 1; i != basis.size(); ++i) {
            if (basis[i].first == basis[i - 1].first) {
                beta[++b] = coasts[basis[i].first][basis[i].second] - alpha[a];
            }
            else if (basis[i].second == basis[i - 1].second) {
                alpha[++a] = coasts[basis[i].first][basis[i].second] - beta[b];
            }
        }
    }

    void sort(vector<pair<size_t, size_t>>& basis) {
        std::sort(basis.begin(), basis.end(), [](pair<size_t, size_t> f, pair<size_t, size_t> s) {
            return (f.first + f.second) < (s.first + s.second);
            });
    }

    vector<vector<T>> sub_matrix(vector<MT> alpha, vector<MT> beta) {
        vector<vector<T>> G(coasts);

        for (size_t i = 1; i != coasts.size(); ++i) {
            for (size_t j = 1; j != coasts[i].size(); ++j) {
                G[i][j] = alpha[i] + beta[j] - coasts[i][j];
            }
        }

        return G;
    }
public:
    void print(vector<vector<T>> v) {
        const int W = 4;
        std::cout.setf(std::cout.right);
        std::cout.fill(' ');
        for (size_t i = 0; i != v.size(); ++i) {
            for (size_t j = 0; j != v[i].size(); ++j) {
                std::cout << std::setw(W) << v[i][j];
            }
            std::cout << '\n';
        }
    }
    void print(vector<pair<size_t, size_t>> v) {
        const int W = 4;
        std::cout.setf(std::cout.right);
        std::cout.fill(' ');
        for (size_t j = 0; j != v.size(); ++j) {
            std::cout << std::setw(W) << v[j].first << std::setw(W) << v[j].second << '\n';
        }
    }
    void print(vector<MT> v) {
        const int W = 4;
        std::cout.setf(std::cout.right);
        std::cout.fill(' ');
        for (size_t j = 0; j != v.size(); ++j) {
            std::cout << std::setw(W) << v[j] << '\n';
        }
    }
    
    T_task(vector<T> producer, vector<T> consumer, vector<vector<T>> coasts) {
        this->coasts = coasts;
        t_body.resize(consumer.size(), vector<T>(producer.size()));
        t_help_body.resize(consumer.size(), vector<T>(producer.size()));
        
        for (size_t i = 1; i != producer.size(); ++i) {
            t_body[0][i] = producer[i];
        }
        for (size_t i = 1; i != consumer.size(); ++i) {
            t_body[i][0] = consumer[i];
        }

        vector<pair<size_t, size_t>> basis = find_base_sln(producer, consumer);
        append_basis(basis);

        vector<MT> alpha(consumer.size(), 0), beta(producer.size(), 0);

        find_potential(alpha, beta, basis);

        vector<vector<T>> G = sub_matrix(alpha, beta);

        std::cout << "\n\nAlpha\n";
        print(alpha);
        std::cout << "\n\nBeta\n";
        print(beta);
        std::cout << "\n\nBody\n";
        print(t_body);
        std::cout << "\n\nBasis\n";
        print(basis);
        std::cout << "\n\nG\n";
        print(G);
        std::cout << "\n\nCoasts\n";
        print(coasts);
    } 
};

int main()
{
    //vector<T> producer{ 0, 5, 3, 5, 1 }; // 4
    //vector<T> consumer{ 0, 3, 2, 1, 4, 2, 2 }; // 6
    
    vector<T> producer{ 0, 70, 30, 20, 40 }; // 4
    vector<T> consumer{ 0, 90, 30, 40 }; // 3

    T_task t(producer, consumer, { 
        {0, 0, 0, 0, 0},
        {0, 2, 3, 4, 3},
        {0, 5, 3, 1, 2},
        {0, 2, 1, 4, 2}
        });
}