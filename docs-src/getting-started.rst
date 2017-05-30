

Getting started
=================

Clone the repository
-----------------------

.. code:: bash

    git clone ...


Setup Docker
----------------

Install Docker
^^^^^^^^^^^^^^^^

.. code:: bash

    ~$ sudo apt-get update
    ~$ sudo apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D
    ~$ sudo apt-get update
    ~$ sudo apt-add-repository 'deb https://apt.dockerproject.org/repo ubuntu-xenial main'
    ~$ sudo apt-get update
    ~$ sudo apt-cache policy docker-engine
    ~$ apt-cache policy docker-engine
    ~$ sudo apt-get install -y docker-engine
    ~$ sudo systemctl status docker

Add current user to docker group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash
	  
   ~$ sudo usermod -aG docker $(whoami)
   ~$ docker ps

You will likely need to restart to ensure that the second command does not throw an error.

Build docker image
^^^^^^^^^^^^^^^^^^^^^

.. code:: bash
	  
  ~$ docker build -t asgaines/arvos .

