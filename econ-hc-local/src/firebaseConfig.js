import firebase from 'firebase'
import 'firebase/firestore'

// Your web app's Firebase configuration
var firebaseConfig = {
	apiKey: "AIzaSyCg749u5HKVAVJ_tJqmnD1PWNNfTeQUBd4",
	authDomain: "econ-hc.firebaseapp.com",
	databaseURL: "https://econ-hc.firebaseio.com",
	projectId: "econ-hc",
	storageBucket: "econ-hc.appspot.com",
	messagingSenderId: "762717017371",
	appId: "1:762717017371:web:ebc9b937dd58e28efc2fbe"
};
// Initialize Firebase
firebase.initializeApp(firebaseConfig);

const db = firebase.firestore();
const auth = firebase.auth();

const usersCollection = db.collection('users')
const pairsCollection = db.collection('pairs')

export {
	db,
	auth,
	usersCollection,
	pairsCollection
}