#!/usr/bin/env python
import os
import matplotlib.pyplot as plt
from automan.api import Problem, Automator


class OneCase(Problem):
    def get_name(self):
        return 'test_case_1_two_spheres_normal_impact'

    def get_commands(self):
        return [
            ('glass', 'python 1_test_elastic_normal_impact_of_two_identical_spheres/two_spheres_normal.py', None),
        ]

    def run(self):
        self.make_output_dir()


class TwoCase(Problem):
    def get_name(self):
        return 'test_case_2_sphere_with_plane_normal_impact'

    def get_commands(self):
        return [
            ('al_alloy', 'python 2_test_elastic_normal_impact_of_sphere_with_plane/sphere_plane_normal_impact.py', None),
        ]

    def run(self):
        self.make_output_dir()

if __name__ == '__main__':
    automator = Automator(
        simulation_dir='outputs',
        output_dir='figures',
        all_problems=[OneCase, TwoCase]
    )
    automator.run()
