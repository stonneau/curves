import os
import unittest

from numpy import array_equal, isclose, matrix, random
from numpy.linalg import norm

from curves import (bezier, bezier_from_hermite, bezier_from_polynomial, cubic_hermite_spline, curve_constraints,
                    exact_cubic, hermite_from_bezier, hermite_from_polynomial, piecewise_bezier_curve,
                    piecewise_cubic_hermite_curve, piecewise_polynomial_curve, polynomial, polynomial_from_bezier,
                    polynomial_from_hermite)


class TestCurves(unittest.TestCase):
    # def print_str(self, inStr):
    #   print inStr
    #   return

    def test_bezier(self):
        print("test_bezier")
        # To test :
        # - Functions : constructor, min, max, derivate,compute_derivate, compute_primitive
        # - Variables : degree, nbWayPoints
        __EPS = 1e-6
        waypoints = matrix([[1., 2., 3.]]).T
        a = bezier(waypoints, 0., 2.)
        t = 0.
        while t < 2.:
            self.assertTrue(norm(a(t) - matrix([1., 2., 3.]).T) < __EPS)
            t += 0.1
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        # time_waypoints = matrix([0., 1.]).transpose()
        # Create bezier6 and bezier
        a = bezier(waypoints, 0., 3.)
        # Test waypoints
        self.assertTrue(a.nbWaypoints == 2)
        for i in range(0, a.nbWaypoints):
            if i == 0:
                self.assertTrue((a.waypointAtIndex(0) == matrix([1., 2., 3.]).transpose()).all())
            elif i == 1:
                self.assertTrue((a.waypointAtIndex(1) == matrix([4., 5., 6.]).transpose()).all())
        # self.assertTrue((a.waypoints == waypoints).all())
        # Test : Degree, min, max, derivate
        # self.print_str(("test 1")
        self.assertEqual(a.degree, a.nbWaypoints - 1)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a(a.min()) == matrix([1., 2., 3.]).transpose()).all())
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        a = a.compute_derivate(100)
        prim = a.compute_primitive(1)
        # Check primitive and derivate - order 1
        for i in range(10):
            t = float(i) / 10.
            self.assertTrue((a(t) == prim.derivate(t, 1)).all())
        self.assertTrue((prim(0) == matrix([0., 0., 0.])).all())
        # Check primitive and derivate - order 2
        prim = a.compute_primitive(2)
        for i in range(10):
            t = float(i) / 10.
            self.assertTrue((a(t) == prim.derivate(t, 2)).all())
        self.assertTrue((prim(0) == matrix([0., 0., 0.])).all())
        # Create new bezier curve
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.], [4., 5., 6.]]).transpose()
        a0 = bezier(waypoints)
        a1 = bezier(waypoints, 0., 3.)
        prim0 = a0.compute_primitive(1)
        prim1 = a1.compute_primitive(1)
        # Check change in argument time_t of bezier
        for i in range(10):
            t = float(i) / 10.
            self.assertTrue(norm(a0(t) - a1(3 * t)) < __EPS)
            self.assertTrue(norm(a0.derivate(t, 1) - a1.derivate(3 * t, 1) * 3.) < __EPS)
            self.assertTrue(norm(a0.derivate(t, 2) - a1.derivate(3 * t, 2) * 9.) < __EPS)
            self.assertTrue(norm(prim0(t) - prim1(t * 3) / 3.) < __EPS)
        self.assertTrue((prim(0) == matrix([0., 0., 0.])).all())
        # testing bezier with constraints
        c = curve_constraints()
        c.init_vel = matrix([0., 1., 1.]).transpose()
        c.end_vel = matrix([0., 1., 1.]).transpose()
        c.init_acc = matrix([0., 1., -1.]).transpose()
        c.end_acc = matrix([0., 100., 1.]).transpose()
        # Check derivate with constraints
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = bezier(waypoints, c)
        self.assertTrue(norm(a.derivate(0, 1) - c.init_vel) < 1e-10)
        self.assertTrue(norm(a.derivate(1, 2) - c.end_acc) < 1e-10)

        # Test serialization : bezier 3
        a.saveAsText("serialization_curve.test")
        # waypoints = matrix([[0,0,0,], [0,0,0,]]).transpose()
        b = bezier()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")

        # Bezier dim 4
        waypoints = matrix([[1., 2., 3., 4.]]).T
        a = bezier(waypoints, 0., 2.)
        # Test serialization : bezier of dim 4
        a.saveAsText("serialization_curve.test")
        # waypoints = matrix([[0,0,0,], [0,0,0,]]).transpose()
        b = bezier()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        return

    def test_polynomial(self):
        print("test_polynomial")
        # To test :
        # - Functions : constructor, min, max, derivate, serialize, deserialize
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = polynomial(waypoints)  # Defined on [0.,1.]
        a = polynomial(waypoints, -1., 3.)  # Defined on [-1.,3.]
        a.min()
        a.max()
        a(0.4)
        # Test get coefficient at degree
        self.assertTrue((a.coeff() == waypoints).all())
        self.assertTrue((a.coeffAtDegree(0) == matrix([1., 2., 3.]).transpose()).all())
        self.assertTrue((a.coeffAtDegree(1) == matrix([4., 5., 6.]).transpose()).all())
        # Other tests
        self.assertTrue((a(a.min()) == matrix([1., 2., 3.]).transpose()).all())
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        a_derivated = a.compute_derivate(1)
        self.assertTrue((a.derivate(0.4, 1) == a_derivated(0.4)).all())
        # Test serialization
        a.saveAsText("serialization_curve.test")
        b = polynomial()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        return

    def test_polynomial_from_boundary_condition(self):
        p0 = matrix([1., 3., -2.]).T
        p1 = matrix([0.6, 2., 2.5]).T
        dp0 = matrix([-6., 2., -1.]).T
        dp1 = matrix([10., 10., 10.]).T
        ddp0 = matrix([1., -7., 4.5]).T
        ddp1 = matrix([6., -1., -4]).T
        min = 1.
        max = 2.5
        polC0 = polynomial(p0, p1, min, max)
        self.assertEqual(polC0.min(), min)
        self.assertEqual(polC0.max(), max)
        self.assertTrue(array_equal(polC0(min), p0))
        self.assertTrue(array_equal(polC0(max), p1))
        self.assertTrue(array_equal(polC0((min + max) / 2.), 0.5 * p0 + 0.5 * p1))
        polC1 = polynomial(p0, dp0, p1, dp1, min, max)
        self.assertEqual(polC1.min(), min)
        self.assertEqual(polC1.max(), max)
        self.assertTrue(isclose(polC1(min), p0).all())
        self.assertTrue(isclose(polC1(max), p1).all())
        self.assertTrue(isclose(polC1.derivate(min, 1), dp0).all())
        self.assertTrue(isclose(polC1.derivate(max, 1), dp1).all())
        polC2 = polynomial(p0, dp0, ddp0, p1, dp1, ddp1, min, max)
        self.assertEqual(polC2.min(), min)
        self.assertEqual(polC2.max(), max)
        self.assertTrue(isclose(polC2(min), p0).all())
        self.assertTrue(isclose(polC2(max), p1).all())
        self.assertTrue(isclose(polC2.derivate(min, 1), dp0).all())
        self.assertTrue(isclose(polC2.derivate(max, 1), dp1).all())
        self.assertTrue(isclose(polC2.derivate(min, 2), ddp0).all())
        self.assertTrue(isclose(polC2.derivate(max, 2), ddp1).all())
        # check that the exception are correctly raised :
        try:
            polC0 = polynomial(p0, p1, max, min)
            self.assertTrue(False)  # should never get there
        except ValueError:
            pass

        try:
            polC1 = polynomial(p0, dp0, p1, dp1, max, min)
            self.assertTrue(False)  # should never get there
        except ValueError:
            pass

        try:
            polC2 = polynomial(p0, dp0, ddp0, p1, dp1, ddp1, max, min)
            self.assertTrue(False)  # should never get there
        except ValueError:
            pass

        return

    def test_cubic_hermite_spline(self):
        print("test_cubic_hermite_spline")
        points = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        tangents = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        time_points = matrix([0., 1.]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        # Test serialization
        a.saveAsText("serialization_curve.test")
        b = cubic_hermite_spline()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        # test dim 4
        points = matrix([[1., 2., 3., 4.], [4., 5., 6., 7.]]).transpose()
        tangents = matrix([[1., 2., 3., 4.], [4., 5., 6., 7.]]).transpose()
        time_points = matrix([0., 1.]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        # Test serialization
        a.saveAsText("serialization_curve.test")
        b = cubic_hermite_spline()
        b.loadFromText("serialization_curve.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_curve.test")
        return

    def test_piecewise_polynomial_curve(self):
        print("test_piecewise_polynomial_curve")
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous, serialize, deserialize
        waypoints0 = matrix([[0., 0., 0.]]).transpose()
        waypoints1 = matrix([[1., 1., 1.]]).transpose()
        waypoints2 = matrix([[1., 1., 1.], [1., 1., 1.]]).transpose()
        polynomial(waypoints0, 0., 0.1)
        a = polynomial(waypoints1, 0., 1.)
        b = polynomial(waypoints2, 1., 3.)
        pc = piecewise_polynomial_curve(a)
        pc.append(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc(pc.min()) == matrix([1., 1., 1.]).transpose()).all())
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise_polynomial_curve()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_piecewise_from_points_list(self):
        N = 7
        points = matrix(random.rand(3, N))
        points_derivative = matrix(random.rand(3, N))
        points_second_derivative = matrix(random.rand(3, N))
        time_points = matrix(random.rand(N)).T
        time_points.sort(0)
        polC0 = piecewise_polynomial_curve.FromPointsList(points, time_points)
        self.assertEqual(polC0.min(), time_points[0, 0])
        self.assertEqual(polC0.max(), time_points[-1, 0])
        self.assertTrue(polC0.is_continuous(0))
        self.assertTrue(not polC0.is_continuous(1))
        for i in range(N):
            self.assertTrue(isclose(polC0(time_points[i, 0]), points[:, i]).all())

        polC1 = piecewise_polynomial_curve.FromPointsList(points, points_derivative, time_points)
        self.assertEqual(polC1.min(), time_points[0, 0])
        self.assertEqual(polC1.max(), time_points[-1, 0])
        self.assertTrue(polC1.is_continuous(0))
        self.assertTrue(polC1.is_continuous(1))
        self.assertTrue(not polC1.is_continuous(2))
        for i in range(N):
            self.assertTrue(isclose(polC1(time_points[i, 0]), points[:, i]).all())
            self.assertTrue(isclose(polC1.derivate(time_points[i, 0], 1), points_derivative[:, i]).all())

        polC2 = piecewise_polynomial_curve.FromPointsList(points, points_derivative, points_second_derivative,
                                                          time_points)
        self.assertEqual(polC2.min(), time_points[0, 0])
        self.assertEqual(polC2.max(), time_points[-1, 0])
        self.assertTrue(polC2.is_continuous(0))
        self.assertTrue(polC2.is_continuous(1))
        self.assertTrue(polC2.is_continuous(2))
        self.assertTrue(not polC2.is_continuous(3))
        for i in range(N):
            self.assertTrue(isclose(polC2(time_points[i, 0]), points[:, i]).all())
            self.assertTrue(isclose(polC2.derivate(time_points[i, 0], 1), points_derivative[:, i]).all())
            self.assertTrue(isclose(polC2.derivate(time_points[i, 0], 2), points_second_derivative[:, i]).all())

        # check if exepetion are corectly raised when time_points are not in ascending values
        time_points[0, 0] = 1
        time_points[1, 0] = 0.5
        try:
            polC0 = piecewise_polynomial_curve.FromPointsList(points, time_points)
            self.assertTrue(False)  # should not get here
        except ValueError:
            pass
        try:
            polC1 = piecewise_polynomial_curve.FromPointsList(points, points_derivative, time_points)
            self.assertTrue(False)  # should not get here
        except ValueError:
            pass
        try:
            polC2 = piecewise_polynomial_curve.FromPointsList(points, points_derivative, points_second_derivative,
                                                              time_points)
            self.assertTrue(False)  # should not get here
        except ValueError:
            pass
        return

    def test_piecewise_bezier3_curve(self):
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = bezier(waypoints, 0., 1.)
        b = bezier(waypoints, 1., 2.)
        pc = piecewise_bezier_curve(a)
        pc.add_curve(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc(pc.min()) == matrix([1., 2., 3.]).transpose()).all())
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise_bezier_curve()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_piecewise_cubic_hermite_curve(self):
        print("test_piecewise_cubic_hermite_curve")
        # To test :
        # - Functions : constructor, min, max, derivate, add_curve, is_continuous
        points = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        tangents = matrix([[2., 2., 2.], [4., 4., 4.]]).transpose()
        time_points0 = matrix([0., 1.]).transpose()
        time_points1 = matrix([1., 2.]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points0)
        b = cubic_hermite_spline(points, tangents, time_points1)
        pc = piecewise_cubic_hermite_curve(a)
        pc.add_curve(b)
        pc.min()
        pc.max()
        pc(0.4)
        self.assertTrue((pc(0.) == matrix([1., 2., 3.]).transpose()).all())
        self.assertTrue((pc.derivate(0.4, 0) == pc(0.4)).all())
        pc.derivate(0.4, 2)
        pc.is_continuous(0)
        pc.is_continuous(1)
        # Test serialization
        pc.saveAsText("serialization_pc.test")
        pc_test = piecewise_cubic_hermite_curve()
        pc_test.loadFromText("serialization_pc.test")
        self.assertTrue((pc(0.4) == pc_test(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_exact_cubic(self):
        print("test_exact_cubic")
        # To test :
        # - Functions : constructor, min, max, derivate, getNumberSplines, getSplineAt
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        time_waypoints = matrix([0., 1.]).transpose()
        a = exact_cubic(waypoints, time_waypoints)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        a.getNumberSplines()
        a.getSplineAt(0)
        # Test serialization
        a.saveAsText("serialization_pc.test")
        b = exact_cubic()
        b.loadFromText("serialization_pc.test")
        self.assertTrue((a(0.4) == b(0.4)).all())
        os.remove("serialization_pc.test")
        return

    def test_exact_cubic_constraint(self):
        print("test_exact_cubic_constraint")
        # To test :
        # - Functions : constructor, min, max, derivate
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        time_waypoints = matrix([0., 1.]).transpose()
        c = curve_constraints()
        c.init_vel = matrix([0., 1., 1.]).transpose()
        c.end_vel = matrix([0., 1., 1.]).transpose()
        c.init_acc = matrix([0., 1., 1.]).transpose()
        c.end_acc = matrix([0., 1., 1.]).transpose()
        c.init_vel
        c.end_vel
        c.init_acc
        c.end_acc
        exact_cubic(waypoints, time_waypoints)
        exact_cubic(waypoints, time_waypoints, c)
        return

    def test_cubic_hermite_spline_2(self):
        points = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        tangents = matrix([[2., 2., 2.], [4., 4., 4.]]).transpose()
        time_points = matrix([0., 1.]).transpose()
        a = cubic_hermite_spline(points, tangents, time_points)
        a.min()
        a.max()
        a(0.4)
        self.assertTrue((a(0.) == matrix([1., 2., 3.]).transpose()).all())
        self.assertTrue((a.derivate(0., 1) == matrix([2., 2., 2.]).transpose()).all())
        self.assertTrue((a.derivate(0.4, 0) == a(0.4)).all())
        a.derivate(0.4, 2)
        return

    def test_conversion_curves(self):
        print("test_conversion_curves")
        __EPS = 1e-6
        waypoints = matrix([[1., 2., 3.], [4., 5., 6.]]).transpose()
        a = bezier(waypoints)
        # converting bezier to polynomial
        a_pol = polynomial_from_bezier(a)
        self.assertTrue(norm(a(0.3) - a_pol(0.3)) < __EPS)
        # converting polynomial to hermite
        a_chs = hermite_from_polynomial(a_pol)
        self.assertTrue(norm(a_chs(0.3) - a_pol(0.3)) < __EPS)
        # converting hermite to bezier
        a_bc = bezier_from_hermite(a_chs)
        self.assertTrue(norm(a_chs(0.3) - a_bc(0.3)) < __EPS)
        self.assertTrue(norm(a(0.3) - a_bc(0.3)) < __EPS)
        # converting bezier to hermite
        a_chs = hermite_from_bezier(a)
        self.assertTrue(norm(a(0.3) - a_chs(0.3)) < __EPS)
        # converting hermite to polynomial
        a_pol = polynomial_from_hermite(a_chs)
        self.assertTrue(norm(a_pol(0.3) - a_chs(0.3)) < __EPS)
        # converting polynomial to bezier
        a_bc = bezier_from_polynomial(a_pol)
        self.assertTrue(norm(a_bc(0.3) - a_pol(0.3)) < __EPS)
        self.assertTrue(norm(a(0.3) - a_bc(0.3)) < __EPS)
        return


if __name__ == '__main__':
    unittest.main()
