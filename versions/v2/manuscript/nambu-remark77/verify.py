import argparse
import hashlib
import json
from itertools import combinations, product
from pathlib import Path

import sympy as sp


def cremona(vector):
    return sp.Matrix([sp.prod(vector[column] for column in range(4) if column != row)
                      for row in range(4)])


def lift(matrix, vector):
    return sp.Matrix.vstack(matrix.adjugate() * vector[:3, :], sp.Matrix([matrix.det()]))


def check_zero(label, expressions, records):
    residuals = [sp.expand(expression) for expression in expressions]
    if any(residual != 0 for residual in residuals):
        raise AssertionError(label)
    records.append({'identity': label, 'zero_residuals': len(residuals)})


def verify_siblings(records):
    eigenvalues = sp.symbols('mu1:4')
    scale = sp.symbols('c')
    variables = sp.symbols('y1:4')
    last = sp.symbols('w')
    first_form = sum(variable ** 2 for variable in variables)
    second_form = sum(eigenvalues[index] * variables[index] ** 2 for index in range(3))
    third_form = last ** 2 + scale * sum(
        sp.prod(eigenvalues[other] for other in range(3) if other != index)
        * variables[index] ** 2 for index in range(3))
    members = []
    for index in range(3):
        others = [other for other in range(3) if other != index]
        expected_first = sum((eigenvalues[other] - eigenvalues[index])
                             * variables[other] ** 2 for other in others)
        members.append(-eigenvalues[index] * first_form + second_form - expected_first)
        expected_second = last ** 2 + scale * sp.prod(
            eigenvalues[index] - eigenvalues[other] for other in others) * variables[index] ** 2
        members.append(-scale * eigenvalues[index] * sum(eigenvalues[other] for other in others)
                       * first_form + scale * eigenvalues[index] * second_form
                       + third_form - expected_second)
    check_zero('six reducible net members', members, records)
    difference = [eigenvalues[1] - eigenvalues[2], eigenvalues[2] - eigenvalues[0],
                  eigenvalues[0] - eigenvalues[1]]
    vandermonde = sp.prod(eigenvalues[row] - eigenvalues[column]
                          for row in range(3) for column in range(row + 1, 3))
    net_coefficients = sp.Matrix([
        difference + [0],
        [eigenvalues[index] * difference[index] for index in range(3)] + [0],
        [scale * sp.prod(eigenvalues[other] for other in range(3) if other != index)
         * difference[index] for index in range(3)] + [-scale * vandermonde]])
    check_zero('normalized diagonal net coefficient sums',
               net_coefficients * sp.ones(4, 1), records)
    check_zero('equal-square scaled base coordinates',
               [sp.prod(difference) + vandermonde], records)
    points = [tuple(signs) + (1,) for signs in product([-1, 1], repeat=3)]
    planes = [(first, second, sign) for first, second in combinations(range(4), 2)
              for sign in [-1, 1]]
    incidence = [[point[first] == sign * point[second] for first, second, sign in planes]
                 for point in points]
    assert all(sum(row) == 6 for row in incidence)
    assert all(sum(row[index] for row in incidence) == 4 for index in range(12))
    hadamard = sp.Matrix([[-1, -1, 1, 1], [1, -1, -1, 1],
                          [-1, 1, -1, 1], [1, 1, 1, 1]])
    normal = 2 * sp.eye(4) - sp.ones(4)

    def projective_key(vector):
        pivot = next(entry for entry in vector if entry != 0)
        return tuple(sp.cancel(entry / pivot) for entry in vector)

    observed = {projective_key(hadamard * sp.Matrix(point)) for point in points}
    expected = {projective_key(sp.eye(4).col(index)) for index in range(4)}
    expected |= {projective_key(normal.col(index)) for index in range(4)}
    assert observed == expected and len(observed) == 8
    spectral, step = sp.symbols('a h')
    first, second, third = variables
    coordinates = sp.Matrix([first, second, third, last])
    field = sp.Matrix([0, spectral * first * third, -spectral * first * second])
    collision = lift(last * sp.eye(3) - step * field.jacobian(variables) / 2, coordinates)
    transformed = sp.Matrix([collision[0], collision[1] + sp.I * collision[2],
                              collision[1] - sp.I * collision[2], collision[3]])
    plus = last + sp.I * step * spectral * first / 2
    minus = last - sp.I * step * spectral * first / 2
    expected_collision = sp.Matrix([first * plus * minus,
                                    (second + sp.I * third) * minus ** 2,
                                    (second - sp.I * third) * plus ** 2,
                                    last * plus * minus])
    check_zero('double-spectrum cubic lift', transformed - expected_collision, records)
    positive, negative = sp.symbols('zplus zminus')
    diagonal_coordinates = sp.Matrix([first, positive, negative, last])

    def iterate(index):
        return sp.Matrix([first * (plus * minus) ** index, positive * minus ** (2 * index),
                          negative * plus ** (2 * index), last * (plus * minus) ** index])

    one_step = iterate(1)
    degrees = []
    for index in range(4):
        current = iterate(index)
        composed = one_step.subs(dict(zip(diagonal_coordinates, current)), simultaneous=True)
        check_zero(f'double-spectrum composition n={index}',
                   composed - (plus * minus) ** (2 * index) * iterate(index + 1), records)
        polynomials = [sp.Poly(expression, *diagonal_coordinates, extension=sp.I)
                       for expression in iterate(index + 1).subs({spectral: 2, step: 1})]
        common = polynomials[0]
        for polynomial in polynomials[1:]:
            common = sp.gcd(common, polynomial)
        assert common.total_degree() == 0
        degree = max(polynomial.total_degree() for polynomial in polynomials)
        assert degree == 2 * (index + 1) + 1
        degrees.append(degree)
    scalars = sp.symbols('p q r s')
    left, right = sp.Matrix(sp.symbols('u0:3')), sp.Matrix(sp.symbols('v0:3'))
    check_zero('Hamiltonian pencil reparametrization',
               (scalars[0] * left + scalars[1] * right).cross(
                   scalars[2] * left + scalars[3] * right)
               - (scalars[0] * scalars[3] - scalars[1] * scalars[2]) * left.cross(right), records)
    return {'sign_points': len(points), 'planes': len(planes),
            'incidences': sum(sum(row) for row in incidence), 'planes_per_point': 6,
            'points_per_plane': 4, 'reducible_quadric_members': 6,
            'normal_form_base_point_images': len(observed),
            'double_spectrum_degrees_n_1_through_4': degrees,
            'finite_scope': 'All sign points and planes are enumerated; four collision iterates '
                            'are checked explicitly. The general degree laws and the exhaustive '
                            'six-member classification have written proofs in note.tex.'}


def verify():
    records = []
    alpha, beta, gamma, step = sp.symbols('alpha beta gamma h')
    first, second, third, last = sp.symbols('x y z w')
    coordinates = sp.Matrix([first, second, third, last])
    coefficients = [alpha ** 2, beta ** 2, gamma ** 2]
    euler_field = sp.Matrix([coefficients[0] * second * third,
                             coefficients[1] * first * third,
                             coefficients[2] * first * second])
    jacobian = euler_field.jacobian([first, second, third])
    euler_matrix = last * sp.eye(3) - step * jacobian / 2
    euler_lift = lift(euler_matrix, coordinates)
    transform = sp.Matrix([
        [-beta * gamma * step, -alpha * gamma * step, alpha * beta * step, -2],
        [beta * gamma * step, -alpha * gamma * step, -alpha * beta * step, -2],
        [-beta * gamma * step, alpha * gamma * step, -alpha * beta * step, -2],
        [beta * gamma * step, alpha * gamma * step, alpha * beta * step, -2]])
    normal = 2 * sp.eye(4) - sp.ones(4)
    check_zero('Euler cubic factorization',
               transform * euler_lift + normal * cremona(transform * coordinates) / 8,
               records)
    check_zero('determinant of T',
               [transform.det() - 32 * alpha ** 2 * beta ** 2 * gamma ** 2 * step ** 3],
               records)
    check_zero('normal involution', list(normal ** 2 - 4 * sp.eye(4)), records)
    check_zero('determinant of normal matrix', [normal.det() + 16], records)
    check_zero('Cremona composition',
               cremona(cremona(coordinates)) - sp.prod(coordinates) ** 2 * coordinates,
               records)
    first_scaled, second_scaled, third_scaled, fourth_scaled = sp.symbols('r s t v')
    scaled = sp.Matrix([first_scaled, second_scaled, third_scaled, fourth_scaled])
    hadamard = sp.Matrix([[-1, -1, 1, 1], [1, -1, -1, 1],
                          [-1, 1, -1, 1], [1, 1, 1, 1]])
    quartic = sum(entry ** 4 for entry in scaled)
    quartic -= 2 * sum(scaled[row] ** 2 * scaled[column] ** 2
                       for row in range(4) for column in range(row + 1, 4))
    quartic += 8 * sp.prod(scaled)
    check_zero('Hadamard quartic product', [sp.prod(hadamard * scaled) - quartic], records)
    coordinate_scale = sp.diag(step * beta * gamma, step * alpha * gamma,
                               step * alpha * beta, -2)
    scaled_values = coordinate_scale * coordinates
    cubics = []
    for row in range(3):
        others = [column for column in range(3) if column != row]
        cubics.append(scaled_values[row]
                      * (scaled_values[3] ** 2 + sum(scaled_values[column] ** 2 for column in others)
                         - scaled_values[row] ** 2)
                      - 2 * sp.prod(scaled_values[column] for column in others) * scaled_values[3])
    cubics.append(scaled_values[3] ** 3
                  - scaled_values[3] * sum(scaled_values[column] ** 2 for column in range(3))
                  + 2 * sp.prod(scaled_values[column] for column in range(3)))
    check_zero('four scaled adjugate cubics',
               coordinate_scale * euler_lift - sp.Matrix(cubics) / 4, records)
    matrix = sp.Matrix(3, 3, sp.symbols('p0:9'))
    left = sp.Matrix(sp.symbols('u0:3'))
    right = sp.Matrix(sp.symbols('v0:3'))
    check_zero('cross-product covariance',
               (matrix * left).cross(matrix * right) - matrix.adjugate().T * left.cross(right),
               records)
    eigenvalues = sp.symbols('mu1:4')
    difference_vector = [eigenvalues[1] - eigenvalues[2],
                         eigenvalues[2] - eigenvalues[0],
                         eigenvalues[0] - eigenvalues[1]]
    vandermonde = sp.prod(eigenvalues[row] - eigenvalues[column]
                          for row in range(3) for column in range(row + 1, 3))
    weighted = sum(sp.prod(eigenvalues[column] for column in range(3) if column != row)
                   * difference_vector[row] for row in range(3))
    check_zero('eight-point net identities',
               [sum(difference_vector),
                sum(eigenvalues[row] * difference_vector[row] for row in range(3)),
                weighted - vandermonde], records)
    negative_transform = transform.subs(step, -step)
    check_zero('step reversal matrix identity',
               list(negative_transform + normal * transform / 2), records)
    frame = sp.Matrix([[2, 1, 0], [0, 1, 1], [0, 0, 1]])
    inverse = frame.inv()
    first_form = inverse.T * inverse
    second_form = inverse.T * sp.diag(1, 2, 5) * inverse
    nambu_field = 4 * (first_form * coordinates[:3, :]).cross(second_form * coordinates[:3, :])
    nambu_jacobian = nambu_field.jacobian([first, second, third])
    nambu_lift = lift(last * sp.eye(3) - sp.Rational(3, 10) * nambu_jacobian, coordinates)
    homogeneous_frame = sp.diag(frame, 1)
    special_transform = transform.subs({alpha: sp.sqrt(6), beta: 2 * sp.sqrt(2) * sp.I,
                                        gamma: sp.sqrt(2), step: sp.Rational(3, 5)})
    conjugator = special_transform * homogeneous_frame.inv()
    check_zero('nonorthogonal positive-definite example',
               conjugator * nambu_lift + normal * cremona(conjugator * coordinates) / 8,
               records)
    check_zero('adjugate congruence in the nonorthogonal example',
               list(frame.T * first_form * second_form.adjugate() * first_form * frame
                    - sp.diag(10, 5, 2) / frame.det() ** 2), records)
    negative_lift = lift(last * sp.eye(3) + sp.Rational(3, 10) * nambu_jacobian, coordinates)
    inverse_composition = negative_lift.subs(dict(zip(coordinates, nambu_lift)), simultaneous=True)
    forward_composition = nambu_lift.subs(dict(zip(coordinates, negative_lift)), simultaneous=True)
    check_zero('inverse composition constant in the nonorthogonal example',
               inverse_composition - sp.prod(conjugator * coordinates) ** 2 * coordinates / 256,
               records)
    check_zero('forward composition constant in the nonorthogonal example',
               forward_composition - sp.prod(normal * conjugator * coordinates) ** 2
               * coordinates / 65536, records)
    if first_form * second_form == second_form * first_form:
        raise AssertionError('Example must have noncommuting quadratic forms')
    for leading_size in range(1, 4):
        if first_form[:leading_size, :leading_size].det() <= 0:
            raise AssertionError('First form must be positive definite')
        if second_form[:leading_size, :leading_size].det() <= 0:
            raise AssertionError('Second form must be positive definite')
    siblings = verify_siblings(records)
    return {'schema_version': 1, 'arithmetic': 'exact symbolic polynomial arithmetic',
            'sympy_version': sp.__version__, 'status': 'all identities equal zero',
            'identities': records,
            'sibling_checks': siblings,
            'example': {'P': [list(frame.row(row)) for row in range(3)],
                        'generalized_eigenvalues': [1, 2, 5], 'h': '3/5',
                        'positive_definite_forms': True, 'commuting_forms': False},
            'scope': 'Universal polynomial identities and one exact nonorthogonal example; '
                     'spectral reduction, the general complex-pencil and collision results, '
                     'projective consequences and distinctness arguments '
                     'are proved in note.tex.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    result = verify()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(result, indent=2, default=str) + '\n'
    if arguments.output:
        arguments.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
