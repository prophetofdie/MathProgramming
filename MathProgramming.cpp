#define NOMINMAX
#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <iomanip>
#include <Windows.h>
#include <string>
#include <limits>

using std::vector;
using std::array;
enum class EDirection { Column, Row };

using T = double;
template<size_t rows, size_t cols>
using matrix = std::array<std::array<T, cols>, rows>;

using node = std::pair<size_t, size_t>;


constexpr std::streamsize alignment = 4;

enum colors { green = 2, blue = 11, red = 12, yellow = 14 };

void set_color(T num, colors color, std::ostream& out) {
    // выбор цвета (в значение color подставлять из перечисления colors)
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
    out << num;
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 15);
}

template<typename T, size_t r, size_t c>
void print(const array<array<T, c>, r>& v) {
    std::cout.setf(std::cout.right);
    std::cout.fill(' ');
    for (auto i = v.begin(); i != v.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j) {
            std::cout << std::setw(alignment) << *j;
        }
        std::cout << '\n';
    }
}
template<size_t x>
void print(const array<node, x>& v) {
    std::cout.setf(std::cout.right);
    std::cout.fill(' ');
    for (auto j = v.begin(); j != v.end(); ++j) {
        std::cout << std::setw(alignment) << j->first << std::setw(alignment) << j->second << '\n';
    }
}
template<size_t x>
void print(const array<T, x>& v) {
    std::cout.setf(std::cout.right);
    std::cout.fill(' ');
    for (auto j = v.begin(); j != v.end(); ++j) {
        std::cout << std::setw(alignment) << *j << '\n';
    }
}
template<size_t r, size_t c, size_t x>
void print_color(const matrix<r, c>& v, const array<node, x>& chain) {
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
template<size_t r, size_t c>
consteval array<node, r + c - 1> find_base_sln(matrix<r, c>& t_body, const array<T, c>& producer, array<T, r> consumer) {
    array<node, r + c - 1> basis_elems;
    T ostatok = producer[0];
    for (size_t i = 0, j = 0, k = 0; ;) {
        if (ostatok - consumer[i] > 0) {
            t_body[i][j] = consumer[i];
            basis_elems[k++] = { i, j };
            ostatok -= consumer[i];
            ++i;
        }
        else if (ostatok - consumer[i] < 0) {
            t_body[i][j] = ostatok;
            basis_elems[k++] = { i, j };
            consumer[i] -= ostatok;
            ++j;
            ostatok = producer[j];
        }
        else {//ostatok == 0
            t_body[i][j] = ostatok;
            basis_elems[k++] = { i, j };
            ++j;
            ++i;
            if (j >= producer.size()) break;
            ostatok = producer[j];
        }
    }
    return basis_elems;
}

template<typename Ty>
consteval Ty abs_(Ty a) {
    return (a < 0) ? -a : a;
}

template<size_t x, size_t y>
consteval array<T, x> gauss_elimination(matrix<x, y> coefficients) {
    size_t n = x;

    // Прямой ход Гаусса
    for (size_t i = 0; i != n; ++i) {
        // Находим строку с максимальным элементом в текущем столбце
        size_t maxRow = i;
        for (size_t k = i + 1; k != n; ++k) {
            if (abs_(coefficients[k][i]) > abs_(coefficients[maxRow][i])) {
                maxRow = k;
            }
        }

        // Меняем местами строки
        std::swap(coefficients[maxRow], coefficients[i]);

        // Приводим главную диагональ к единице
        for (size_t j = i + 1; j != n; ++j) {
            T ratio = coefficients[j][i] / coefficients[i][i];
            for (size_t k = i; k != n + 1; ++k) {
                coefficients[j][k] -= coefficients[i][k] * ratio;
            }
        }
    }

    // Обратный ход для нахождения решения
    array<T, x> solution;
    for (size_t i = n; i != 0; --i) {
        solution[i - 1] = coefficients[i - 1][n] / coefficients[i - 1][i - 1];
        for (size_t k = i - 1; k != 0; --k) {
            coefficients[k - 1][n] -= coefficients[k - 1][i - 1] * solution[i - 1];
        }
    }

    return solution;  // Возвращаем вектор решений
}

template<size_t r, size_t c>
consteval void find_potential(const matrix<r, c>& coasts, array<T, r>& alpha, array<T, c>& beta, const array<node, r + c - 1>& basis) {
    const size_t basis_size = basis.size();
    const size_t coast_size = coasts.size();
    const size_t alpha_size = alpha.size();

    matrix<r + c, r + c + 1> A;
    array<T, r + c + 1> d;

    d.fill(0);
    A.fill(d);
    for (size_t i = 0; i != basis_size; ++i) {
        d[i] = coasts[basis[i].first][basis[i].second];
    }

    for (size_t i = 0; i != basis_size; ++i) {

        A[i][basis[i].first] = true;
        A[i][basis[i].second + coast_size] = true;
        A[i].back() = d[i];
    }

    d.back() = 0;// обнуление Alpha
    A.back()[basis[0].first] = true;// обнуление Alpha

    //std::cout << "\n\nA matrix\n";
    //print(A);

    auto solution = gauss_elimination(A);

    size_t i = 0;
    for (size_t k = 0; i != alpha_size; ++i, ++k) {
        alpha[k] = solution[i];
    }
    for (size_t k = 0; i != solution.size(); ++i, ++k) {
        beta[k] = solution[i];
    }
}

template<size_t r, size_t c>
consteval matrix<r, c> find_G_matrix(const matrix<r, c>& coasts, const array<T, r>& alpha, const array<T, c>& beta) {
    matrix<r, c> G(coasts);

    for (size_t i = 0; i != r; ++i) {
        for (size_t j = 0; j != c; ++j) {
            G[i][j] = alpha[i] + beta[j] - coasts[i][j];
        }
    }

    return G;
}

template<size_t r, size_t c, size_t ch>
consteval T find_min_in_chain(matrix<r, c>& t_body, const array<node, ch>& chain) {
    T min = std::numeric_limits<T>::max();
    for (auto i = chain.begin(); i != chain.end(); i++) {
        if (min > t_body[i->first][i->second] and t_body[i->first][i->second] != 0) {
            min = t_body[i->first][i->second];
        }
    }
    return min;
}

template<size_t x, size_t ch>
consteval bool chain_recurs(const array<node, x>& basis, array<node, ch>& chain, EDirection direction, size_t i = 0) {
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

    ++i;
    if (i != 1 and get_axis(chain[i - 1]) == get_axis(chain[0])) {
        return true;
    }
    else {
        for (decltype(auto) n : basis) {
            if (n != chain[i - 1] and get_axis(chain[i - 1]) == get_axis(n)) {
                chain[i] = n;
                if (chain_recurs(basis, chain, (direction == EDirection::Column) ? EDirection::Row : EDirection::Column, i)) {
                    return true;
                }
                else {
                    --i;
                }
            }
        }
    }
    return false;
}

template<size_t x, size_t y>
consteval array<node, y> compress_array__(const array<node, x>& a) {
    array<node, y> res;
    for (size_t i = 0; i != y; ++i) {
        res[i] = a[i];
    }
    return res;
}

template<size_t x>
consteval auto compress_array(const array<node, x>& a) {
    size_t size = 0;
    // Loop through the array to find the index of the sentinel value
    for (size_t i = 0; i < x; ++i) {
        if (a[i] == node({ -1, -1 })) {
            break;
        }
        ++size;
    }
    // Pass the constant size to the next function
    return compress_array__<x, 4>(a);
}

/// <summary>
/// Составление цепочки
/// </summary>
/// <param name="basis">Кооринаты элементов базиса</param>
/// <param name="destination">Первый элемент цепочки</param>
/// <returns></returns>
template<size_t x>
consteval auto find_chain(const array<node, x>& basis, const node& destination) {
    array<node, x + 1> chain;
    chain.fill({ -1,-1 });
    chain[0] = destination;

    if (chain_recurs(basis, chain, EDirection::Row)) {
        return compress_array(chain);
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
template<size_t r, size_t c>
consteval node find_new_place(const matrix<r, c>& G) {
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
template<size_t r, size_t c, size_t x>
consteval node balance(matrix<r, c>& t_body, const array<node, x>& chain, const node& new_basis_place) {
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
template<size_t r, size_t c>
consteval bool is_optimal(const matrix<r, c>& G) {
    for (auto i = G.begin(); i != G.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j) {
            if (*j > 0) return false;
        }
    }
    return true;
}

template<size_t r, size_t c>
consteval array<node, r + c - 1> base_sln(const array<T, c>& producer, const array<T, r>& consumer, matrix<r, c>& t_body) {
    array<node, r + c - 1> basis = find_base_sln(t_body, producer, consumer);

    //Дополняем базис
    for (auto i = basis.begin(); i != basis.end(); ++i) {
        if (not((i->first == (i + 1)->first) or (i->second == (i + 1)->second))) {
            basis.back() = { i->first, i->second + 1 };
            break;
        }
    }
    return basis;
}
template<size_t r, size_t c>
consteval void optimise(matrix<r, c>& t_body, array<node, r + c - 1>& basis, const matrix<r, c>& coasts) {
    //std::cout << "\n\n////////////////////////////////////////////\n";
    //std::cout << "Iteration: " << ++counter << '\n';
    array<T, c> beta;
    array<T, r> alpha;

    find_potential(coasts, alpha, beta, basis);

    matrix G = find_G_matrix(coasts, alpha, beta);

    //std::cout << "\n\nBasis\n";
    //print_color(t_body, basis);
    //std::cout << "\n\nAlpha\n";
    //print(alpha);
    //std::cout << "\n\nBeta\n";
    //print(beta);
    //std::cout << "\n\nG\n";
    //print(G);

    if (is_optimal(G)) return;

    node new_basis_place = find_new_place(G);

    auto chain = find_chain(basis, new_basis_place);

    node aborted = balance(t_body, chain, new_basis_place);

    //std::cout << "\n\nBalanced Body\n";
    //print_color(t_body, chain);

    *std::find(basis.begin(), basis.end(), aborted) = new_basis_place;

    optimise(t_body, basis, coasts);
}

template<size_t r, size_t c>
consteval T T_task(const array<T, c>& producer, const array<T, r>& consumer, const matrix<r, c>& coast) {
    matrix<r, c> t_body;

    {
        array<T, c> temp;
        temp.fill(0);
        t_body.fill(temp);
    }

    array<node, r + c - 1> basis = base_sln(producer, consumer, t_body);

    optimise(t_body, basis, coast);

    T z = 0;
    const size_t i_size = t_body.size(), j_size = t_body[0].size();
    for (size_t i = 0; i != i_size; ++i) {
        for (size_t j = 0; j != j_size; ++j) {
            z += t_body[i][j] * coast[i][j];
        }
    }

    //std::cout << "\n\nOptimal function value: " << z;
    return z;
}

int main() {
    set_color(0, colors::green, std::cout); // Для красоты

    constexpr size_t r = 3;
    constexpr size_t c = 4;

    constexpr array<T, c> producer{ 70, 30, 20, 40 };
    constexpr array<T, r> consumer{ 90, 30, 40 };

    constexpr matrix<r, c> coast{
        array<T, c>{ 2, 3, 4, 3},
        array<T, c>{ 5, 3, 1, 2},
        array<T, c>{ 2, 1, 4, 2}
    };

    constexpr T z = T_task(producer, consumer, coast);
}