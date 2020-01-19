sudo apt-get -y update
sudo apt-get -y dist-upgrade

sudo apt-get -y install python3-pip
sudo apt-get -y install nginx

sudo pip3 install boto3
sudo pip3 install setuptools
sudo pip3 install waitress
sudo pip3 install flask

sudo openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout /etc/nginx/cert.key -out /etc/nginx/cert.crt
#sudo openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout cert.key -out cert.crt
sudo openssl dhparam -out /etc/nginx/dhparam.pem 4096
#openssl dhparam -out dhparam.pen 4096

sudo cat <<EOF >/etc/nginx/nginx.conf


events { }

http {
    server {
        listen 80;
        server_name ai-real.org www.ai-real.org;

        return 301 https://$host$request_uri;
	    }

	server {
		listen 443;
		server_name ai-real.org www.ai-real.org;

		ssl                  on;
		ssl_certificate      /etc/letsencrypt/live/www.ai-real.org/fullchain.pem;
		ssl_certificate_key  /etc/letsencrypt/live/www.ai-real.org/privkey.pem;

		ssl_prefer_server_ciphers on;
		ssl_session_timeout 1d;
		ssl_session_cache shared:SSL:50m;
		ssl_session_tickets off;

		ssl_dhparam /etc/nginx/dhparam.pem;

		# https://wiki.mozilla.org/Security/Server_Side_TLS
		# Intermediate compatibility (recommended)
		ssl_protocols TLSv1.2 TLSv1.3;
		ssl_ciphers 'TLS_AES_128_GCM_SHA256:TLS_AES_256_GCM_SHA384:TLS_CHACHA20_POLY1305_SHA256:ECDHE-ECDSA-AES128-GCM-SHA256:ECDHE-RSA-AES128-GCM-SHA256:ECDHE-ECDSA-AES256-GCM-SHA384:ECDHE-RSA-AES256-GCM-SHA384:ECDHE-ECDSA-CHACHA20-POLY1305:ECDHE-RSA-CHACHA20-POLY1305:DHE-RSA-AES128-GCM-SHA256:DHE-RSA-AES256-GCM-SHA384';

		location / {
		    proxy_pass http://localhost:8080;
		    proxy_set_header Host $host;
		    proxy_set_header X-Real-IP $remote_addr;
		    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
		    proxy_set_header X-Forwarded-Proto https;
			}
		}
	}

EOF




# CERTBOT
# Requires the following prmissions.
#route53:ListHostedZones
#route53:GetChange
#route53:ChangeResourceRecordSets

sudo apt-get update
sudo apt-get install software-properties-common
sudo add-apt-repository universe
sudo add-apt-repository ppa:certbot/certbot
sudo apt-get update
sudo apt-get install certbot python-certbot-nginx
sudo apt-get install python3-certbot-dns-route53
sudo certbot certonly --dns-route53 -d www.ai-real.org -d ai-real.org -i nginx
(crontab -l 2>/dev/null; echo "*0 12 * * * /usr/bin/certbot renew --quiet") | crontab -



# Install postgresql
sudo apt-get -y install curl ca-certificates gnupg
curl https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -
sudo add-apt-repository 'deb http://apt.postgresql.org/pub/repos/apt/ bionic-pgdg main'
sudo apt-get -y install postgresql-12 postgresql-client-12 postgresql-server-dev-12
sudo pip3 install psycopg2

sudo systemctl stop postgresql@12-main

VARLIB='/var/lib/postgresql/12/main'
ETC='/etc/postgresql/12/main'
sudo -u postgres mv $VARLIB /var/lib/postgresql/12/_main

sudo -u postgres /usr/lib/postgresql/12/bin/initdb -D $VARLIB --wal-segsize=1
sudo -u postgres sed -i "s|max_wal_size = 1GB|max_wal_size = 64MB|" $ETC/postgresql.conf
sudo -u postgres sed -i "s|min_wal_size = 80MB|min_wal_size = 5MB|" $ETC/postgresql.conf

sudo systemctl enable postgresql@12-main
sudo systemctl start postgresql@12-main

sudo -u postgres psql -c 'create database production;'
sudo -u postgres psql -c "create user webserver with encrypted password 'ecdrvd98u';"
sudo -u postgres psql -c 'grant all privileges on database production to webserver;'








