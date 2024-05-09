#pragma once

#include "Collider.h"
#include "omp.h"

namespace Pivot {
class Magnetic {
  private:
    friend class Simulation;

  public:
    void Solve(SurfaceMesh &mesh) {
        m_Mesh = &mesh;

        InitSolver();
        // SolveMagneticByMC();
        SolveMagneticByFPI();
    }

  private:
    struct Sampler {
        std::default_random_engine RandEngine;
        std::discrete_distribution<> SegmentSampler;
        std::uniform_real_distribution<double> RandRussianRoulette;
    };
    void InitSolver() {
        m_MagneticPressure.resize(m_Mesh->size(), 0);
        m_MagneticHn.resize(m_Mesh->size(), 0);
        m_MagneticHt.resize(m_Mesh->size(), 0);
    }
    void SolveMagneticByFPI() {
        int size = m_Mesh->size();
        VectorXd u(size);
        VectorXd utmp(size);
        VectorXd b(size);
        MatrixXd A(size, size);

        for (int i = 0; i < size; i++) {
            b(i) = -2 * m_Lambda * m_Hext.dot(m_Mesh->Normals[i]);
            u(i) = b(i) / (1 - m_Lambda);
            for (int j = 0; j < size; j++) {
                if (i == j) {
                    A(i, j) = 0;
                    continue;
                }
                // A(i, j) =
                //     2 * m_Lambda *
                //     dGdxdClamped(m_Mesh->Positions[i], m_Mesh->Positions[j],
                //                  m_Mesh->Normals[i], m_EpsFPI) *
                //     m_Mesh->Areas[j];
                A(i, j) = 2 * m_Lambda *
                          dGdxd(m_Mesh->Positions[i], m_Mesh->Positions[j],
                                m_Mesh->Normals[i], m_EpsFPI * m_EpsFPI) *
                          m_Mesh->Areas[j];
            }
        }

        // fmt::print("\n");
        int iter;
        for (iter = 0; iter < m_NumIteration; iter++) {
            utmp = A * u + b;
            double L1 = (u - utmp).cwiseAbs().sum() / size;
            double maxCoeff = (u - utmp).cwiseAbs().maxCoeff();
            // fmt::print("Iter [{:02d}] L1({:.5f}) maxCoeff({:.5f})\n", iter,
            // L1,
            //            maxCoeff);
            u = utmp;
            if (maxCoeff < m_StopThres) {
                break;
            }
        }
        utmp = A * u + b;
        double L1 = (u - utmp).cwiseAbs().sum() / size;
        double maxCoeff = (u - utmp).cwiseAbs().maxCoeff();
        // fmt::print("\nIter [{:02d}] L1({:.5e}) maxCoeff({:.5e})\n", iter, L1,
        //            maxCoeff);
        for (int i = 0; i < size; i++) {
            double w = m_MU * (1 + m_Chi) / (-m_Chi) * u(i);
            m_MagneticHn[i] = -w / (m_MU * (1 + m_Chi));
            double Hn_ = m_MagneticHn[i] * (1 + m_Chi);

            Vector2d nx = m_Mesh->Normals[i];
            Vector2d tx = Vector2d(nx.y(), -nx.x());
            m_MagneticHt[i] = m_Hext.dot(tx);
            for (int j = 0; j < size; j++) {
                if (i == j) {
                    continue;
                }
                // m_MagneticHt[i] +=
                //     -dGdxdClamped(m_Mesh->Positions[i], m_Mesh->Positions[j],
                //                   tx, m_EpsFPI) *
                //     m_Mesh->Areas[j] * u(j);
                m_MagneticHt[i] +=
                    -dGdxd(m_Mesh->Positions[i], m_Mesh->Positions[j], tx,
                           m_EpsFPI * m_EpsFPI) *
                    m_Mesh->Areas[j] * u(j);
            }

            double pressure = 0;

            pressure += m_MU * (1 + m_Chi) * 0.5 *
                        (m_MagneticHn[i] * m_MagneticHn[i] -
                         m_MagneticHt[i] * m_MagneticHt[i]);
            pressure -=
                m_MU * 0.5 * (Hn_ * Hn_ - m_MagneticHt[i] * m_MagneticHt[i]);
            m_MagneticPressure[i] = pressure;
        }
    }
    void SolveMagneticByMC() {
        SetRandomEngine();
        int size = m_MagneticPressure.size();
#pragma omp parallel for schedule(dynamic) default(shared)
        for (int i = 0; i < size; i++) {
            int tid = omp_get_thread_num();
            int xIdx = i;
            Vector2d x = m_Mesh->Positions[xIdx];
            Vector2d nx = m_Mesh->Normals[xIdx];
            Vector2d tx = Vector2d(nx.y(), -nx.x());
            double sumHn = 0;
            double sumHt = 0;
            for (int sample = 0; sample < m_NumSample; sample++) {
                double weight;
                double value;
                int yIdx = _UniformSample(xIdx, tid);
                Vector2d y = m_Mesh->Positions[yIdx];
                value = WoB(yIdx, tid);

                weight = dGdxd(x, y, nx, m_EpsMC) * m_Mesh->TotalArea;
                sumHn += -weight * value - m_Lambda * weight * value +
                         m_Lambda * m_Hext.dot(nx);
                sumHn += m_Hext.dot(nx);

                weight = dGdxd(x, y, tx, m_EpsMC) * m_Mesh->TotalArea;
                sumHt += -weight * value;
                sumHt += m_Hext.dot(tx);
            }
            m_MagneticHn[i] = sumHn / m_NumSample;
            m_MagneticHt[i] = sumHt / m_NumSample;

            double pressure = 0;

            pressure += m_MU * (1 + m_Chi) * 0.5 *
                        (m_MagneticHn[i] * m_MagneticHn[i] -
                         m_MagneticHt[i] * m_MagneticHt[i]);
            double Hn_ = m_MagneticHn[i] * (1 + m_Chi);
            pressure -=
                m_MU * 0.5 * (Hn_ * Hn_ - m_MagneticHt[i] * m_MagneticHt[i]);
            m_MagneticPressure[i] = pressure;
        }
    }
    void SetRandomEngine() {
        std::random_device rdevice;
        int max_num_threads = omp_get_max_threads();
        m_Samplers.resize(max_num_threads);
        for (int i = 0; i < max_num_threads; i++) {
            m_Samplers[i].RandEngine = std::default_random_engine(rdevice());
            m_Samplers[i].SegmentSampler = std::discrete_distribution<>(
                m_Mesh->Areas.begin(), m_Mesh->Areas.end());
            m_Samplers[i].RandRussianRoulette =
                std::uniform_real_distribution<double>(0, 1);
        }
    }
    double WoB(int idx, int tid = 0) {
        Vector2d x = m_Mesh->Positions[idx];
        Vector2d nx = m_Mesh->Normals[idx];

        double boundaryValue = 2 * m_Lambda * m_Hext.dot(nx);
        double rr =
            m_Samplers[tid].RandRussianRoulette(m_Samplers[tid].RandEngine);
        if (rr > m_RussianRoulette) {
            return -boundaryValue;
        }

        int yIdx = _UniformSample(idx, tid);
        Vector2d y = m_Mesh->Positions[yIdx];
        double weight = 2 * m_Lambda * dGdxd(x, y, nx, m_EpsMC) *
                        m_Mesh->TotalArea / m_RussianRoulette;
        return weight * WoB(yIdx, tid) - boundaryValue;
    }
    int _UniformSample(int xIdx, int tid = 0) {
        int idx;
        do {
            idx = m_Samplers[tid].SegmentSampler(m_Samplers[tid].RandEngine);
        } while (idx == xIdx);
        return idx;
    }

    double dGdxd(const Vector2d &x, const Vector2d &y, const Vector2d &xd,
                 double eps) {
        Vector2d r = y - x;
        return 1.0 / (2.0 * m_PI) * r.dot(xd) /
               (std::max)(r.squaredNorm(), eps);
    }
    double dGdxdClamped(const Vector2d &x, const Vector2d &y,
                        const Vector2d &xd, double eps) {
        Vector2d r = y - x;
        return (r.norm() < eps)
                   ? 0
                   : (1.0 / (2.0 * m_PI) * r.dot(xd) / r.squaredNorm());
    }

  private:
    inline static const double m_PI = 3.141592653589783;
    inline static const double m_MU = 4e-7 * m_PI;

    SurfaceMesh *m_Mesh;
    std::vector<double> m_MagneticPressure;
    std::vector<double> m_MagneticHn;
    std::vector<double> m_MagneticHt;

    double m_Chi = .5;
    double m_Lambda = (-m_Chi) / (2 + m_Chi);
    Vector2d m_Hext = Vector2d(0, 5e4);

    double m_RussianRoulette = 0.5;
    int m_NumSample = 20000;
    double m_EpsMC = 1e-6;
    std::vector<Sampler> m_Samplers;

    int m_NumIteration = 20;
    double m_EpsFPI = 1e-3;
    double m_StopThres = 1e-6;
};
} // namespace Pivot
