import numpy as np

"""Grade average for Msc in Physics and BSc in Engineering Physics
"""

master = {
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

bachelor = {
    'Linear algebra and geometry': {
        'hp': 4.5,
        'grade': 3,
    },
    'Introductory mathematical analysis': {
        'hp': 6,
        'grade': 4,
    },
    'Real analysis': {
        'hp': 6,
        'grade': 4,
    },
    'Computer programming': {
        'hp': 6,
        'grade': 4,
    },
    'Multivariable analysis': {
        'hp': 6,
        'grade': 5,
    },
    'Mechanics 1': {
        'hp': 7.5,
        'grade': 5,
    },
    'Linear algebra and numerical analysis': {
        'hp': 7.5,
        'grade': 5,
    },
    'Mechanics 2': {
        'hp': 6,
        'grade': 4,
    },
    'Vector fields and classical physics': {
        'hp': 4.5,
        'grade': 4,
    },
    'Complex mathematical analysis': {
        'hp': 6,
        'grade': 4,
    },
    'Electrical circuits and systems': {
        'hp': 7.5,
        'grade': 5,
    },
    'Theory of electromagnetic fields': {
        'hp': 7.5,
        'grade': 5,
    },
    'Optics': {
        'hp': 6.0,
        'grade': 5,
    },
    'Fourier analysis': {
        'hp': 6.0,
        'grade': 5,
    },
    'Mathematical statistics': {
        'hp': 4.5,
        'grade': 5,
    },
    'Automatic control': {
        'hp': 4.5,
        'grade': 5,
    },
    'Environmental physics': {
        'hp': 4.5,
        'grade': 5,
    },
    'Experimental physics 1 - measuring technique': {
        'hp': 9,
        'grade': 5,
    },
    'Quantum physics': {
        'hp': 6.0,
        'grade': 5,
    },
    'Thermodynamics and statistical mechanics': {
        'hp': 7.5,
        'grade': 5,
    },
    'Applied quantum physics': {
        'hp': 4.5,
        'grade': 5,
    },
    'Special relativity': {
        'hp': 4.5,
        'grade': 5,
    },
    'Fluid mechanics': {
        'hp': 4.5,
        'grade': 5,
    },
    'Solid state physics': {
        'hp': 7.5,
        'grade': 5,
    },
    'Subatomic physics': {
        'hp': 6.0,
        'grade': 5,
    },
    'Experimental physics 2 - basis': {
        'hp': 4.5,
        'grade': 5,
    },
    "Bachelor's thesis in Physics": {
        'hp': 15,
        'grade': 5,
    },
}

courses = {**bachelor, **master}

cumWeightedGrade = 0
cumHp = 0
for course, data in courses.items():
    hp = data['hp']
    grade = data['grade']
    cumWeightedGrade += hp*grade
    cumHp += hp
    print(course, data['grade'])
weightedGrade = cumWeightedGrade / cumHp

print(f'Weighted grade average: {weightedGrade:.2f}')