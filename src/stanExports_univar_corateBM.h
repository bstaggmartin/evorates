// Generated by rstantools.  Do not edit by hand.

/*
    backwards-BM-simulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    backwards-BM-simulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with backwards-BM-simulator.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_univar_corateBM_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_univar_corateBM");
    reader.add_event(71, 69, "end", "model_univar_corateBM");
    return reader;
}
#include <stan_meta_header.hpp>
class model_univar_corateBM : public prob_grad {
private:
        int n;
        int n_e;
        vector_d X;
        matrix_d eV;
        vector_d prune_T;
        std::vector<std::vector<int> > des_e;
        std::vector<int> tip_e;
        std::vector<int> real_e;
        std::vector<int> prune_seq;
        double R0_prior;
        double Rsig2_prior;
        double X0_prior;
        matrix_d chol_eV;
public:
    model_univar_corateBM(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_univar_corateBM(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_univar_corateBM_namespace::model_univar_corateBM";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "n_e", "int", context__.to_vec());
            n_e = int(0);
            vals_i__ = context__.vals_i("n_e");
            pos__ = 0;
            n_e = vals_i__[pos__++];
            current_statement_begin__ = 5;
            validate_non_negative_index("X", "n", n);
            context__.validate_dims("data initialization", "X", "vector_d", context__.to_vec(n));
            X = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                X(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 6;
            validate_non_negative_index("eV", "n_e", n_e);
            validate_non_negative_index("eV", "n_e", n_e);
            context__.validate_dims("data initialization", "eV", "matrix_d", context__.to_vec(n_e,n_e));
            eV = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n_e, n_e);
            vals_r__ = context__.vals_r("eV");
            pos__ = 0;
            size_t eV_j_2_max__ = n_e;
            size_t eV_j_1_max__ = n_e;
            for (size_t j_2__ = 0; j_2__ < eV_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < eV_j_1_max__; ++j_1__) {
                    eV(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("prune_T", "((2 * n) - 1)", ((2 * n) - 1));
            context__.validate_dims("data initialization", "prune_T", "vector_d", context__.to_vec(((2 * n) - 1)));
            prune_T = Eigen::Matrix<double, Eigen::Dynamic, 1>(((2 * n) - 1));
            vals_r__ = context__.vals_r("prune_T");
            pos__ = 0;
            size_t prune_T_j_1_max__ = ((2 * n) - 1);
            for (size_t j_1__ = 0; j_1__ < prune_T_j_1_max__; ++j_1__) {
                prune_T(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("des_e", "((2 * n) - 1)", ((2 * n) - 1));
            validate_non_negative_index("des_e", "2", 2);
            context__.validate_dims("data initialization", "des_e", "int", context__.to_vec(((2 * n) - 1),2));
            des_e = std::vector<std::vector<int> >(((2 * n) - 1), std::vector<int>(2, int(0)));
            vals_i__ = context__.vals_i("des_e");
            pos__ = 0;
            size_t des_e_k_0_max__ = ((2 * n) - 1);
            size_t des_e_k_1_max__ = 2;
            for (size_t k_1__ = 0; k_1__ < des_e_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < des_e_k_0_max__; ++k_0__) {
                    des_e[k_0__][k_1__] = vals_i__[pos__++];
                }
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("tip_e", "n", n);
            context__.validate_dims("data initialization", "tip_e", "int", context__.to_vec(n));
            tip_e = std::vector<int>(n, int(0));
            vals_i__ = context__.vals_i("tip_e");
            pos__ = 0;
            size_t tip_e_k_0_max__ = n;
            for (size_t k_0__ = 0; k_0__ < tip_e_k_0_max__; ++k_0__) {
                tip_e[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 12;
            validate_non_negative_index("real_e", "n_e", n_e);
            context__.validate_dims("data initialization", "real_e", "int", context__.to_vec(n_e));
            real_e = std::vector<int>(n_e, int(0));
            vals_i__ = context__.vals_i("real_e");
            pos__ = 0;
            size_t real_e_k_0_max__ = n_e;
            for (size_t k_0__ = 0; k_0__ < real_e_k_0_max__; ++k_0__) {
                real_e[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 13;
            validate_non_negative_index("prune_seq", "(n - 1)", (n - 1));
            context__.validate_dims("data initialization", "prune_seq", "int", context__.to_vec((n - 1)));
            prune_seq = std::vector<int>((n - 1), int(0));
            vals_i__ = context__.vals_i("prune_seq");
            pos__ = 0;
            size_t prune_seq_k_0_max__ = (n - 1);
            for (size_t k_0__ = 0; k_0__ < prune_seq_k_0_max__; ++k_0__) {
                prune_seq[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 16;
            context__.validate_dims("data initialization", "R0_prior", "double", context__.to_vec());
            R0_prior = double(0);
            vals_r__ = context__.vals_r("R0_prior");
            pos__ = 0;
            R0_prior = vals_r__[pos__++];
            current_statement_begin__ = 17;
            context__.validate_dims("data initialization", "Rsig2_prior", "double", context__.to_vec());
            Rsig2_prior = double(0);
            vals_r__ = context__.vals_r("Rsig2_prior");
            pos__ = 0;
            Rsig2_prior = vals_r__[pos__++];
            current_statement_begin__ = 18;
            context__.validate_dims("data initialization", "X0_prior", "double", context__.to_vec());
            X0_prior = double(0);
            vals_r__ = context__.vals_r("X0_prior");
            pos__ = 0;
            X0_prior = vals_r__[pos__++];
            // initialize transformed data variables
            current_statement_begin__ = 23;
            validate_non_negative_index("chol_eV", "n_e", n_e);
            validate_non_negative_index("chol_eV", "n_e", n_e);
            chol_eV = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n_e, n_e);
            stan::math::fill(chol_eV, DUMMY_VAR__);
            // execute transformed data statements
            current_statement_begin__ = 24;
            stan::math::assign(chol_eV, cholesky_decompose(eV));
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 28;
            num_params_r__ += 1;
            current_statement_begin__ = 29;
            num_params_r__ += 1;
            current_statement_begin__ = 30;
            num_params_r__ += 1;
            current_statement_begin__ = 31;
            validate_non_negative_index("raw_R", "n_e", n_e);
            num_params_r__ += n_e;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_univar_corateBM() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 28;
        if (!(context__.contains_r("R0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable R0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("R0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "R0", "double", context__.to_vec());
        double R0(0);
        R0 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(R0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable R0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 29;
        if (!(context__.contains_r("Rsig2")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable Rsig2 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("Rsig2");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "Rsig2", "double", context__.to_vec());
        double Rsig2(0);
        Rsig2 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, Rsig2);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable Rsig2: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 30;
        if (!(context__.contains_r("X0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable X0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("X0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "X0", "double", context__.to_vec());
        double X0(0);
        X0 = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(X0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable X0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 31;
        if (!(context__.contains_r("raw_R")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable raw_R missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("raw_R");
        pos__ = 0U;
        validate_non_negative_index("raw_R", "n_e", n_e);
        context__.validate_dims("parameter initialization", "raw_R", "vector_d", context__.to_vec(n_e));
        Eigen::Matrix<double, Eigen::Dynamic, 1> raw_R(n_e);
        size_t raw_R_j_1_max__ = n_e;
        for (size_t j_1__ = 0; j_1__ < raw_R_j_1_max__; ++j_1__) {
            raw_R(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(raw_R);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable raw_R: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 28;
            local_scalar_t__ R0;
            (void) R0;  // dummy to suppress unused var warning
            if (jacobian__)
                R0 = in__.scalar_constrain(lp__);
            else
                R0 = in__.scalar_constrain();
            current_statement_begin__ = 29;
            local_scalar_t__ Rsig2;
            (void) Rsig2;  // dummy to suppress unused var warning
            if (jacobian__)
                Rsig2 = in__.scalar_lb_constrain(0, lp__);
            else
                Rsig2 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 30;
            local_scalar_t__ X0;
            (void) X0;  // dummy to suppress unused var warning
            if (jacobian__)
                X0 = in__.scalar_constrain(lp__);
            else
                X0 = in__.scalar_constrain();
            current_statement_begin__ = 31;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> raw_R;
            (void) raw_R;  // dummy to suppress unused var warning
            if (jacobian__)
                raw_R = in__.vector_constrain(n_e, lp__);
            else
                raw_R = in__.vector_constrain(n_e);
            // transformed parameters
            current_statement_begin__ = 35;
            validate_non_negative_index("R", "n_e", n_e);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> R(n_e);
            stan::math::initialize(R, DUMMY_VAR__);
            stan::math::fill(R, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 36;
            stan::math::assign(R, add(R0, multiply(multiply(stan::math::sqrt(Rsig2), chol_eV), raw_R)));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 35;
            size_t R_j_1_max__ = n_e;
            for (size_t j_1__ = 0; j_1__ < R_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(R(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: R" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable R: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            {
            current_statement_begin__ = 40;
            validate_non_negative_index("SS", "((2 * n) - 1)", ((2 * n) - 1));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> SS(((2 * n) - 1));
            stan::math::initialize(SS, DUMMY_VAR__);
            stan::math::fill(SS, DUMMY_VAR__);
            current_statement_begin__ = 41;
            validate_non_negative_index("XX", "((2 * n) - 1)", ((2 * n) - 1));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> XX(((2 * n) - 1));
            stan::math::initialize(XX, DUMMY_VAR__);
            stan::math::fill(XX, DUMMY_VAR__);
            current_statement_begin__ = 42;
            validate_non_negative_index("VV", "((2 * n) - 1)", ((2 * n) - 1));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> VV(((2 * n) - 1));
            stan::math::initialize(VV, DUMMY_VAR__);
            stan::math::fill(VV, DUMMY_VAR__);
            current_statement_begin__ = 43;
            validate_non_negative_index("LL", "(n - 1)", (n - 1));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> LL((n - 1));
            stan::math::initialize(LL, DUMMY_VAR__);
            stan::math::fill(LL, DUMMY_VAR__);
            current_statement_begin__ = 44;
            int counter(0);
            (void) counter;  // dummy to suppress unused var warning
            stan::math::fill(counter, std::numeric_limits<int>::min());
            current_statement_begin__ = 45;
            validate_non_negative_index("des_X", "2", 2);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> des_X(2);
            stan::math::initialize(des_X, DUMMY_VAR__);
            stan::math::fill(des_X, DUMMY_VAR__);
            current_statement_begin__ = 46;
            validate_non_negative_index("des_V", "2", 2);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> des_V(2);
            stan::math::initialize(des_V, DUMMY_VAR__);
            stan::math::fill(des_V, DUMMY_VAR__);
            current_statement_begin__ = 49;
            lp_accum__.add(cauchy_log<propto__>(R0, 0, R0_prior));
            current_statement_begin__ = 50;
            lp_accum__.add(cauchy_log<propto__>(Rsig2, 0, Rsig2_prior));
            current_statement_begin__ = 51;
            lp_accum__.add(normal_log<propto__>(X0, 0, X0_prior));
            current_statement_begin__ = 52;
            lp_accum__.add(std_normal_log<propto__>(raw_R));
            current_statement_begin__ = 55;
            stan::math::assign(SS, rep_vector(0, ((2 * n) - 1)));
            current_statement_begin__ = 56;
            stan::model::assign(SS, 
                        stan::model::cons_list(stan::model::index_multi(real_e), stan::model::nil_index_list()), 
                        elt_multiply(stan::model::rvalue(prune_T, stan::model::cons_list(stan::model::index_multi(real_e), stan::model::nil_index_list()), "prune_T"), stan::math::exp(R)), 
                        "assigning variable SS");
            current_statement_begin__ = 57;
            stan::model::assign(XX, 
                        stan::model::cons_list(stan::model::index_multi(tip_e), stan::model::nil_index_list()), 
                        X, 
                        "assigning variable XX");
            current_statement_begin__ = 58;
            stan::model::assign(VV, 
                        stan::model::cons_list(stan::model::index_multi(tip_e), stan::model::nil_index_list()), 
                        rep_vector(0, n), 
                        "assigning variable VV");
            current_statement_begin__ = 59;
            stan::math::assign(counter, 0);
            current_statement_begin__ = 60;
            for (auto& i : prune_seq) {
                (void) i;  // dummy to suppress unused var warning
                current_statement_begin__ = 61;
                stan::math::assign(des_X, stan::model::rvalue(XX, stan::model::cons_list(stan::model::index_multi(stan::model::rvalue(des_e, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "des_e")), stan::model::nil_index_list()), "XX"));
                current_statement_begin__ = 62;
                stan::math::assign(des_V, add(stan::model::rvalue(VV, stan::model::cons_list(stan::model::index_multi(stan::model::rvalue(des_e, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "des_e")), stan::model::nil_index_list()), "VV"), stan::model::rvalue(SS, stan::model::cons_list(stan::model::index_multi(stan::model::rvalue(des_e, stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "des_e")), stan::model::nil_index_list()), "SS")));
                current_statement_begin__ = 63;
                stan::math::assign(counter, (counter + 1));
                current_statement_begin__ = 64;
                stan::model::assign(LL, 
                            stan::model::cons_list(stan::model::index_uni(counter), stan::model::nil_index_list()), 
                            (-(0.5) * ((stan::math::log((2 * stan::math::pi())) + stan::math::log(sum(des_V))) + (pow((get_base1(des_X, 1, "des_X", 1) - get_base1(des_X, 2, "des_X", 1)), 2) / sum(des_V)))), 
                            "assigning variable LL");
                current_statement_begin__ = 65;
                stan::model::assign(XX, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (((get_base1(des_V, 2, "des_V", 1) / sum(des_V)) * get_base1(des_X, 1, "des_X", 1)) + ((get_base1(des_V, 1, "des_V", 1) / sum(des_V)) * get_base1(des_X, 2, "des_X", 1))), 
                            "assigning variable XX");
                current_statement_begin__ = 66;
                stan::model::assign(VV, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (1 / ((1 / get_base1(des_V, 1, "des_V", 1)) + (1 / get_base1(des_V, 2, "des_V", 1)))), 
                            "assigning variable VV");
            }
            current_statement_begin__ = 68;
            lp_accum__.add((sum(LL) - (0.5 * ((stan::math::log((2 * stan::math::pi())) + stan::math::log(get_base1(VV, 1, "VV", 1))) + (pow((X0 - get_base1(XX, 1, "XX", 1)), 2) / sum(stan::model::rvalue(VV, stan::model::cons_list(stan::model::index_multi(stan::model::rvalue(des_e, stan::model::cons_list(stan::model::index_uni(1), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "des_e")), stan::model::nil_index_list()), "VV")))))));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("R0");
        names__.push_back("Rsig2");
        names__.push_back("X0");
        names__.push_back("raw_R");
        names__.push_back("R");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_e);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(n_e);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_univar_corateBM_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double R0 = in__.scalar_constrain();
        vars__.push_back(R0);
        double Rsig2 = in__.scalar_lb_constrain(0);
        vars__.push_back(Rsig2);
        double X0 = in__.scalar_constrain();
        vars__.push_back(X0);
        Eigen::Matrix<double, Eigen::Dynamic, 1> raw_R = in__.vector_constrain(n_e);
        size_t raw_R_j_1_max__ = n_e;
        for (size_t j_1__ = 0; j_1__ < raw_R_j_1_max__; ++j_1__) {
            vars__.push_back(raw_R(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 35;
            validate_non_negative_index("R", "n_e", n_e);
            Eigen::Matrix<double, Eigen::Dynamic, 1> R(n_e);
            stan::math::initialize(R, DUMMY_VAR__);
            stan::math::fill(R, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 36;
            stan::math::assign(R, add(R0, multiply(multiply(stan::math::sqrt(Rsig2), chol_eV), raw_R)));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t R_j_1_max__ = n_e;
                for (size_t j_1__ = 0; j_1__ < R_j_1_max__; ++j_1__) {
                    vars__.push_back(R(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_univar_corateBM";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "R0";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Rsig2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "X0";
        param_names__.push_back(param_name_stream__.str());
        size_t raw_R_j_1_max__ = n_e;
        for (size_t j_1__ = 0; j_1__ < raw_R_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "raw_R" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t R_j_1_max__ = n_e;
            for (size_t j_1__ = 0; j_1__ < R_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "R" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "R0";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "Rsig2";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "X0";
        param_names__.push_back(param_name_stream__.str());
        size_t raw_R_j_1_max__ = n_e;
        for (size_t j_1__ = 0; j_1__ < raw_R_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "raw_R" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t R_j_1_max__ = n_e;
            for (size_t j_1__ = 0; j_1__ < R_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "R" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_univar_corateBM_namespace::model_univar_corateBM stan_model;
#endif
