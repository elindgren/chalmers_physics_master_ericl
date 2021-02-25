import numpy as np

"""Grade average for Msc in Physics
"""

courses = {
    'Quantum Mechanics': {
        'hp': 4.5,
        'grade': 4,
    },
    'Learning From Data': {
        'hp': 7.5,
        'grade': 5,
    },
    'Experimental Methods In Modern Physics': {
        'hp': 3.0,
        'grade': 5,
    },
    'Computational physics': {
        'hp': 7.5,
        'grade': 5,
    },
    'Statistical physics': {
        'hp': 7.5,
        'grade': 5,
    },
    'Computational materials and molecular physics': {
        'hp': 7.5,
        'grade': 5,
    },
    'Plasma physics with applications': {
        'hp': 7.5,
        'grade': 5,
    },
    'High performance computing': {
        'hp': 7.5,
        'grade': 5,
    },
    'Quantum computing': {
        'hp': 7.5,
        'grade': 5,
    },
    'Technology and society': {
        'hp': 7.5,
        'grade': 5,
    },
    'Advanced simulation and machine learning': {
        'hp': 7.5,
        'grade': 5,
    },
}

cumWeightedGrade = 0
cumHp = 0
for course, data in courses.items():
    hp = data['hp']
    grade = data['grade']
    cumWeightedGrade += hp*grade
    cumHp += hp
    print(data['grade'])
weightedGrade = cumWeightedGrade / cumHp

print(f'Weighted grade average: {weightedGrade:.2f}')