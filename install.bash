# Run to install needed dependencies system-wide
sudo apt-get update
sudo apt-get upgrade

sudo apt-get install python3-pip
sudo apt-get install libpq-dev

sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common \
    postgresql-client

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo apt-key fingerprint 0EBFCD88
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"

sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io

# Install docker-compose
sudo curl -L "https://github.com/docker/compose/releases/download/1.25.5/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# Permissions
sudo groupadd docker
sudo usermod -aG docker $USER

# Create folders to store docker data
mkdir -p $(pwd)/db/postgresql_data

# This should be done manually (sudo is not working)
sudo echo '{"iptables": false}' > /etc/docker/daemon.json
