#pragma once

#include <tuple>    // import std::tie()
#include <utility>  // import std::move

#include "vector2.hpp"

namespace numeric {

    /**
     * @brief matrix2
     *
     */
    template <typename T1, typename T2 = T1> class matrix2 : public vector2<T1, T2> {
      public:
        /**
         * @brief Construct a new matrix2 object
         *
         * @param x
         * @param y
         */
        constexpr matrix2(T1&& x, T2&& y) noexcept : vector2<T1, T2>{std::move(x), std::move(y)} {}

        // /**
        //  * @brief Construct a new matrix2 object
        //  *
        //  * @param x
        //  * @param y
        //  */
        // constexpr matrix2(const T1& x, const T2& y) : vector2<T1, T2>{x, y} {}

        /** @name Arithmetic operators
         *  definie +, -, *, /, +=, -=, *=, /=, etc.
         */
        ///@{

        /**
         * @brief Negate
         *
         * @return matrix2
         */
        constexpr let operator-() const -> matrix2<T1, T2> { return {-this->x(), -this->y()}; }

        /**
         * @brief Add
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return matrix2&
         */
        template <typename U1, typename U2> constexpr let operator+=(const matrix2<U1, U2>& other)
            -> matrix2<T1, T2>& {
            this->_x += other.x();
            this->_y += other.y();
            return *this;
        }

        /**
         * @brief Substract
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return matrix2&
         */
        template <typename U1, typename U2>  //
        constexpr let operator-=(const matrix2<U1, U2>& other) -> matrix2<T1, T2>& {
            this->_x -= other.x();
            this->_y -= other.y();
            return *this;
        }

        /**
         * @brief Multiply
         *
         * @tparam R
         * @param[in] alpha
         * @return matrix2&
         */
        template <typename R> constexpr let operator*=(const R& alpha) -> matrix2<T1, T2>& {
            this->_x *= alpha;
            this->_y *= alpha;
            return *this;
        }

        /**
         * @brief Divide
         *
         * @tparam R
         * @param[in] alpha
         * @return matrix2&
         */
        template <typename R> constexpr let operator/=(const R& alpha) -> matrix2<T1, T2>& {
            this->_x /= alpha;
            this->_y /= alpha;
            return *this;
        }

        /**
         * @brief Add
         *
         * @tparam U1
         * @tparam U2
         * @param[in] x
         * @param[in] y
         * @return matrix2
         */
        template <typename U1, typename U2>  //
        friend constexpr let operator+(matrix2<T1, T2> x, const matrix2<U1, U2>& y)
            -> matrix2<T1, T2> {
            return std::move(x) += y;
        }

        /**
         * @brief Substract
         *
         * @tparam U1
         * @tparam U2
         * @param[in] x
         * @param[in] y
         * @return matrix2
         */
        template <typename U1, typename U2>  //
        friend constexpr let operator-(matrix2<T1, T2> x, const matrix2<U1, U2>& y)
            -> matrix2<T1, T2> {
            return std::move(x) -= y;
        }

        /**
         * @brief Multiply by a scalar
         *
         * @tparam R
         * @param[in] x
         * @param[in] alpha
         * @return matrix2
         */
        template <typename R> friend constexpr let operator*(matrix2<T1, T2> x, const R& alpha)
            -> matrix2<T1, T2> {
            return x *= alpha;
        }

        /**
         * @brief Multiply (by a scalar)
         *
         * @tparam R
         * @param[in] alpha
         * @param[in] x
         * @return matrix2
         */
        template <typename R> friend constexpr let operator*(const R& alpha, matrix2<T1, T2> x)
            -> matrix2<T1, T2> {
            return x *= alpha;
        }

        /**
         * @brief Divide (by a scalar)
         *
         * @tparam R
         * @param[in] x
         * @param[in] alpha
         * @return matrix2
         */
        template <typename R> friend constexpr let operator/(matrix2<T1, T2> x, const R& alpha)
            -> matrix2<T1, T2> {
            return x /= alpha;
        }

        /**
         * @brief
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return constexpr let
         */
        template <typename U1, typename U2>  //
        constexpr let mdot(const vector2<U1, U2>& other) const -> vector2<U1, U2> {
            return {this->_x.dot(other), this->_y.dot(other)};
        }

        /**
         * @brief
         *
         * @tparam U1
         * @tparam U2
         * @param[in] other
         * @return constexpr let
         */
        constexpr let det() const -> f64 {
            return this->x().x() * this->y().y() - this->x().y() * this->y().x();
        }

        ///@}
    };
}  // namespace numeric
