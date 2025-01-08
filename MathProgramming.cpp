#define NOMINMAX
#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <iomanip>
#include <Windows.h>
#include <string>
#include <limits>
#include <array>

using std::vector;
enum class EDirection { Column, Row };

using T = long long;
using matrix = vector<vector<T>>;

using node = std::pair<size_t, size_t>;


constexpr std::streamsize alignment = 4;

enum colors { green = 2, blue = 11, red = 12, yellow = 14 };

void set_color(T num, colors color, std::ostream& out) { 
    // выбор цвета (в значение color подставлять из перечисления colors)
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
    out << num;
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 15);
}

template<typename T> 
void print(const vector<vector<T>>& v) {
    std::cout.setf(std::cout.right);
    std::cout.fill(' ');
    for (auto i = v.begin(); i != v.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j) {
            std::cout << std::setw(alignment) << *j;
        }
        std::cout << '\n';
    }
}
void print(const vector<node>& v) {
    std::cout.setf(std::cout.right);
    std::cout.fill(' ');
    for (auto j = v.begin(); j != v.end(); ++j) {
        std::cout << std::setw(alignment) << j->first << std::setw(alignment) << j->second << '\n';
    }
}
void print(const vector<T>& v) {
    std::cout.setf(std::cout.right);
    std::cout.fill(' ');
    for (auto j = v.begin(); j != v.end(); ++j) {
        std::cout << std::setw(alignment) << *j << '\n';
    }
}
void print_color(const matrix& v, const vector<node>& chain) {
    std::cout.setf(std::cout.right);
    std::cout.fill(' ');
    const size_t size_i = v.size(), size_j = v[0].size();

    for (size_t i = 0; i != size_i; ++i) {
        for (size_t j = 0; j != size_j; ++j) {
            std::cout << std::setw(alignment);
            if (std::find(chain.begin(), chain.end(), std::make_pair(i, j)) == chain.end()) {
                std::cout << v[i][j];
            }
            else {
                set_color(v[i][j], colors::green, std::cout);
            }
        }
        std::cout << '\n';
    }
}

/// <summary>
/// Нахождение базиса алгоритмом северо-западного угла
/// </summary>
/// <param name="t_body">Матрица количества поставок</param>
/// <param name="producer">Производитель</param>
/// <param name="consumer">Потребитель</param>
/// <returns>Координаты базисных элементов</returns>
vector<node> find_base_sln(matrix& t_body, vector<T> producer, vector<T> consumer) {
    vector<node> basis_elems;
    T ostatok = producer[0];
    for (size_t i = 0, j = 0; ;) {
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

template<typename T>
vector<T> gauss_elimination(matrix& coefficients, vector<T>& rightHandSide) {
    size_t n = coefficients.size();

    // Расширенная матрица (коэффициенты + правая часть)
    for (size_t i = 0; i != n; ++i) {
        coefficients[i].push_back(rightHandSide[i]);
    }

    // Прямой ход Гаусса
    for (size_t i = 0; i != n; ++i) {
        // Находим строку с максимальным элементом в текущем столбце
        size_t maxRow = i;
        for (size_t k = i + 1; k != n; ++k) {
            if (std::abs(coefficients[k][i]) > std::abs(coefficients[maxRow][i])) {
                maxRow = k;
            }
        }

        // Меняем местами строки
        std::swap(coefficients[i], coefficients[maxRow]);

        // Приводим главную диагональ к единице
        for (size_t j = i + 1; j != n; ++j) {
            T ratio = coefficients[j][i] / coefficients[i][i];
            for (size_t k = i; k != n + 1; ++k) {
                coefficients[j][k] -= coefficients[i][k] * ratio;
            }
        }
    }

    // Обратный ход для нахождения решения
    vector<T> solution(n);
    for (size_t i = n; i != 0; --i) {
        solution[i - 1] = coefficients[i - 1][n] / coefficients[i - 1][i - 1];
        for (size_t k = i - 1; k != 0; --k) {
            coefficients[k - 1][n] -= coefficients[k - 1][i - 1] * solution[i - 1];
        }
    }

    return solution;  // Возвращаем вектор решений
}

void find_potential(const matrix& coasts, vector<T>& alpha, vector<T>& beta, const vector<node>& basis) {
    const size_t basis_size = basis.size();
    const size_t coast_size = coasts.size();
    const size_t alpha_size = alpha.size();

    matrix A(basis_size + 1, vector<T>(coasts.size() + coasts[0].size(), 0));
    vector<T> d(basis_size + 1);
    

    for (size_t i = 0; i != basis_size; ++i) {
        d[i] = coasts[basis[i].first][basis[i].second];
    }

    for (size_t i = 0; i != basis_size; ++i) {
        A[i][basis[i].first] = true;
        A[i][basis[i].second + coast_size] = true;
    }

    d.back() = 0;// обнуление Alpha
    A.back()[basis[0].first] = true;// обнуление Alpha

    vector<T> solution = gauss_elimination(A, d);
    
    size_t i = 0;
    for (size_t k = 0; i != alpha_size; ++i, ++k) {
        alpha[k] = solution[i];
    }
    for (size_t k = 0; i != solution.size(); ++i, ++k) {
        beta[k] = solution[i];
    }
}

matrix find_G_matrix(const matrix& coasts, const vector<T>& alpha, const vector<T>& beta) {
    matrix G(coasts);
    const size_t i_size = coasts.size(), j_size = coasts[0].size();

    for (size_t i = 0; i != i_size; ++i) {
        for (size_t j = 0; j != j_size; ++j) {
            G[i][j] = alpha[i] + beta[j] - coasts[i][j];
        }
    }

    return G;
}

T find_min_in_chain(matrix& t_body, const vector<node>& chain) {
    T min = std::numeric_limits<T>::max();
    for (auto i = chain.begin(); i != chain.end(); i++) {
        if (min > t_body[i->first][i->second] and t_body[i->first][i->second] != 0) {
            min = t_body[i->first][i->second];
        }
    }
    return min;
}

bool chain_recurs(const vector<node>& basis, vector<node>& chain, EDirection direction) {
    auto get_axis = [direction](node n) {
        switch (direction) {
        case EDirection::Column:
            return n.first;
        case EDirection::Row:
            return n.second;
        default:
            throw;
        }
        };
    // chain[0] - добавляемый элемент

    if (chain.size() != 1 and get_axis(chain.back()) == get_axis(chain[0])) {
        return true;
    }
    else {
        for (decltype(auto) n : basis) {
            if (n != chain.back() and get_axis(chain.back()) == get_axis(n)) {
                chain.push_back(n);
                if (chain_recurs(basis, chain, (direction == EDirection::Column) ? EDirection::Row : EDirection::Column)) {
                    return true;
                }
                else {
                    chain.pop_back();
                }
            }
        }
    }
    return false;
}

/// <summary>
/// Составление цепочки
/// </summary>
/// <param name="basis">Кооринаты элементов базиса</param>
/// <param name="destination">Первый элемент цепочки</param>
/// <returns></returns>
vector<node> find_chain(const vector<node>& basis, const node& destination) {
    vector<node> chain = { destination }; 

    if (chain_recurs(basis, chain, EDirection::Row)) {
        return chain;
    }
    else {
        throw "Цепочка не собралась";
    }
}

/// <summary>
/// Поиск координат для нового элемента базиса
/// </summary>
/// <param name="G"></param>
/// <returns></returns>
node find_new_place(const matrix& G) {
    const size_t i_size = G.size(), j_size = G[0].size();

    for (size_t i = 0; i != i_size; ++i) {
        for (size_t j = 0; j != j_size; ++j) {
            if (G[i][j] > 0) return { i, j };
        }
    }
    return { -1, -1 };
}

/// <summary>
/// Ввод нового элемента в базис
/// Балансировка матрицы и вывод элемента из базиса
/// </summary>
/// <param name="t_body">Матрица количества поставок</param>
/// <param name="chain">Цепочка</param>
/// <param name="new_basis_place">Координаты нового элемента базиса</param>
/// <returns></returns>
node balance(matrix& t_body, const vector<node>& chain, const node& new_basis_place) {
    T min = find_min_in_chain(t_body, chain);// Находим минимальный элемент для вывода из базиса
    bool sign = true; 
    bool is_abort = false; 
    node aborted; // Выведенный элемент
    for (decltype(auto) chain_el : chain) {
        t_body[chain_el.first][chain_el.second] += min * ((sign) ? 1 : -1);
        sign = !sign;
        if (t_body[chain_el.first][chain_el.second] == 0) {
            is_abort = true; // Вывели элемент
            aborted = { chain_el.first, chain_el.second };
        }
    }
    
    // Если элемент не был выведен -> балансируем еще 
    return (!is_abort) ? balance(t_body, chain, new_basis_place) : aborted;
}

/// <summary>
/// Проверяет решение на оптимальность
/// </summary>
/// <param name="G"></param>
/// <returns>true - если решение оптимальное</returns>
bool is_optimal(const matrix& G) {
    for (auto i = G.begin(); i != G.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j) {
            if (*j > 0) return false;
        }
    }
    return true;
}

vector<node> base_sln(vector<T>& producer, vector<T>& consumer, matrix& t_body) {
    vector<node> basis = find_base_sln(t_body, producer, consumer);

    //Дополняем базис
    for (auto i = basis.begin(); i != basis.end(); ++i) {
        if (not((i->first == (i + 1)->first) or (i->second == (i + 1)->second))) {
            basis.emplace_back(i->first, i->second + 1);
            break;
        }
    }
    return basis;
}

void optimise(matrix& t_body, vector<node>& basis, matrix& coasts, vector<T>& producer, vector<T>& consumer) {
    static size_t counter = 0;
    std::cout << "\n\n////////////////////////////////////////////\n";
    std::cout << "Iteration: " << ++counter << '\n';
    vector<T> beta(producer.size(), 0), alpha(consumer.size(), 0);

    find_potential(coasts, alpha, beta, basis);

    matrix G = find_G_matrix(coasts, alpha, beta);
;
    std::cout << "\n\nBasis\n";
    print_color(t_body, basis);
    std::cout << "\n\nAlpha\n";
    print(alpha);
    std::cout << "\n\nBeta\n";
    print(beta);
    std::cout << "\n\nG\n";
    print(G);
    
    if (is_optimal(G)) return;

    node new_basis_place = find_new_place(G);

    vector<node> chain = find_chain(basis, new_basis_place);
    
    std::cout << "\n\nChained Body\n";
    print_color(t_body, chain);
    
    node aborted = balance(t_body, chain, new_basis_place);
    
    std::cout << "\n\nBalanced Body\n";
    print_color(t_body, chain);
    
    basis.push_back(new_basis_place);
    basis.erase(std::remove_if(basis.begin(), basis.end(), [&aborted](const node& i) { return i == aborted; }));

    optimise(t_body, basis, coasts, producer, consumer);
}

void T_task(vector<T>& producer, vector<T>& consumer, matrix& coast) {
    matrix t_body(consumer.size(), vector<T>(producer.size()));
    vector<node> basis = base_sln(producer, consumer, t_body);

    optimise(t_body, basis, coast, producer, consumer);

    T z = 0;
    const size_t i_size = t_body.size(), j_size = t_body[0].size();
    for (size_t i = 0; i != i_size; ++i) {
        for (size_t j = 0; j != j_size; ++j) {
            z += t_body[i][j] * coast[i][j];
        }
    }

    std::cout << "\n\nOptimal function value: " << z;
}

int main() {
    set_color(0, colors::green, std::cout); // Для красоты
    
    vector<T> producer{ 70, 30, 20, 40 };
    vector<T> consumer{ 90, 30, 40 };

    matrix coast{
        { 2, 3, 4, 3},
        { 5, 3, 1, 2},
        { 2, 1, 4, 2}
    };

    T_task(producer, consumer, coast);
}